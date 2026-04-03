require 'date'
require 'logger'
require 'etna/command'

class Magma
  class RetrieveProjectTemplate < Etna::Command
    string_flags << '--file'
    string_flags << '--target-model'

    def execute(project_name, target_model: 'project', file: "#{project_name}_models_#{target_model}_tree.csv")
      unless ::File.exist?(file)
        puts "File #{file} is being prepared from the #{project_name} project."
        puts "Copying models descending from #{target_model}..."
        prepare_template(file, project_name, target_model)
        puts
        puts "Done!  You can start editing the file #{file} now"
      else
        puts "File #{file} already exists!  Please remove or specify a different file name before running again."
      end
    end

    def workflow
      @workflow ||= Etna::Clients::Magma::AddProjectModelsWorkflow.new(magma_client: magma_client)
    end

    def magma_client
      Etna::Clients::LocalMagmaClient.new
    end

    def prepare_template(file, project_name, target_model)
      tf = Tempfile.new

      begin
        File.open(tf.path, 'wb') { |f| workflow.write_models_template_csv(project_name, target_model, io: f) }
        FileUtils.cp(tf.path, file)
      ensure
        tf.close!
      end
    end

    def setup(config)
      super
      Magma.instance.setup_db
      Magma.instance.load_models
      require_relative './server'
    end
  end

  class LoadProject < Etna::Command
    usage '[project_name, path/to/file.json] # Import attributes into database for given project name from JSON file'

    def execute(project_name, file_name)
      file = File.open(file_name)
      file_data = JSON.parse(file.read, symbolize_names: true)

      file_data[:models].each do |model_name, model_json|
        model_name = model_name.to_s
        template = model_json[:template]

        load_model(project_name, model_name, template)

        template[:attributes].each do |attribute_name, attribute|
          load_attribute(project_name, model_name, attribute)
        end
      end
    end

    def setup(config)
      super
      Magma.instance.setup_db
    end

    private

    def load_model(project_name, model_name, template)
      dictionary_json = if template[:dictionary]
        template[:dictionary][:attributes].merge(
          dictionary_model: template[:dictionary][:dictionary_model]
        )
      else
        nil
      end

      Magma.instance.db[:models].insert(
        project_name: project_name,
        model_name: model_name,
        dictionary: Sequel.pg_json_wrap(dictionary_json),
      )
    end

    def load_attribute(project_name, model_name, attribute)
      row = attribute.
        slice(*options).
        merge(
          project_name: project_name,
          model_name: model_name,
          column_name: attribute[:attribute_name],
          type: attribute[:attribute_type],
          validation: Sequel.pg_json_wrap(attribute[:validation]),
        )

      Magma.instance.db[:attributes].insert(row)
    end

    def options
      require_relative './attribute'
      @options ||= Magma::Attribute.options - [:loader] + [:created_at, :updated_at, :attribute_name]
    end
  end

  class Migrate < Etna::Command
    usage "Run database wide migrations"
    string_flags << '--version'

    def execute(version: nil)
      Sequel.extension(:migration)
      db = Magma.instance.db

      if version
        puts "Migrating to version #{version}, globally"
        Sequel::Migrator.run(db, File.join("db", "migrations"), target: version.to_i)
      else
        puts 'Migrating to latest, globally'
        Sequel::Migrator.run(db, File.join("db", "migrations"))
      end
    end

    def setup(config)
      super
      Magma.instance.setup_db
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected magma instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
      Magma.instance.load_models
    end
  end

  class Load < Etna::Command
    usage 'Run data loaders on models for current dataset.'

    def execute(*args)
      loaders = Magma.instance.find_descendents(Magma::Loader)

      if args.empty?
        # List available loaders
        puts 'Available loaders:'
        loaders.each do |loader|
          puts "%30s  %s" % [ loader.loader_name, loader.description ]
        end
        exit
      end

      loader = loaders.find do |l| l.loader_name == args[0] end

      raise "Could not find a loader named #{args[0]}" unless loader

      loader = loader.new
      loader.load(*args[1..-1])
      begin
        loader.dispatch
      rescue Magma::LoadFailed => e
        puts "Load failed with these complaints:"
        puts e.complaints
      end
    end

    def setup(config)
      super
      Magma.instance.load_models
    end
  end

  class Unload < Etna::Command
    usage '<project_name> <model_name> # Dump the dataset of the model into a tsv'

    def execute(project_name, model_name)
      require_relative './payload'
      require_relative './retrieval'
      require_relative './tsv_writer'

      begin
        model = Magma.instance.get_model(project_name, model_name)
        retrieval = Magma::Retrieval.new(model, 'all', 'all', page: 1, page_size: 100_000)
        payload = Magma::Payload.new
        Magma::TSVWriter.new(model, retrieval, payload).write_tsv{ |lines| puts lines }
      rescue Exception => e
        puts "Unload failed:"
        puts e.message
      end
    end

    def setup(config)
      super
      Magma.instance.load_models
      Magma.instance.setup_db
    end
  end

  class MeasureDataRows < Etna::Command
    def setup(config)
      super
      Magma.instance.load_models
      Magma.instance.setup_db
    end

    def execute
      Magma.instance.magma_projects.keys.each do |project_name|
        project = Magma.instance.get_project(project_name)
        project.models.each do |model_name, model|
          tags = {model_name: model_name.to_s, project_name: project_name.to_s}
          Yabeda.magma.data_rows.set(tags, model.count)
        end
      end
    end
  end

  class AutoDisconnect < Etna::Command
    usage '<project_name> <model_name> # Disconnect all records for a model by setting parent to nil'
    boolean_flags << '--auto'

    def execute(project_name, model_name, auto: false)
      require_relative './server'
      require_relative './loader'

      @project_name = project_name
      model = get_model(model_name)
      
      # Validate model has a parent
      parent_attr = find_parent_attribute(model)
      unless parent_attr
        puts "Error: Model #{model_name} does not have a parent attribute and cannot be disconnected."
        return
      end

      if auto
        disconnect_with_hierarchy(model)
      else
        disconnect_model(model)
      end
    end

    def setup(config)
      super
      Magma.instance.load_models
      Magma.instance.setup_db
    end

    private

    def get_model(model_name)
      Magma.instance.get_model(@project_name, model_name)
    end

    def find_parent_attribute(model)
      model.attributes.values.find { |attr| attr.is_a?(Magma::ParentAttribute) }
    end

    def disconnect_model(model)
      model_name = model.model_name
      parent_attr = find_parent_attribute(model)
      parent_attr_name = parent_attr.name

      # Get all connected records (those with a parent)
      connected_records = model.where(Sequel.~(parent_attr.column_name.to_sym => nil)).all
      
      if connected_records.empty?
        puts "No connected records found for model #{model_name}."
        return
      end

      puts "Found #{connected_records.count} connected record(s) for model #{model_name}."
      puts "Disconnecting records..."

      # Build revisions hash for the loader
      revisions = {}
      revisions[model_name.to_sym] = {}

      connected_records.each do |record|
        identifier = record.send(model.identity.name)
        revisions[model_name.to_sym][identifier.to_sym] = {
          parent_attr_name => nil
        }
      end

      # Use the Magma loader to perform the disconnect
      user = Etna::User.new(
        email: 'system@mountetna.org',
        name: 'System',
        perm: 'a:administration'
      )

      loader = Magma::Loader.new(user, @project_name)

      revisions[model_name.to_sym].each do |record_name, revision|
        loader.push_record(model, record_name.to_s, revision)
      end

      begin
        payload = loader.dispatch_record_set
        puts "Successfully disconnected #{connected_records.count} record(s) from model #{model_name}."
      rescue Magma::LoadFailed => e
        puts "Disconnect failed with errors:"
        puts e.complaints
      rescue Exception => e
        puts "Disconnect failed:"
        puts e.message
      end
    end

    def disconnect_with_hierarchy(model)
      # Use similar pattern to ReparentModelAction to walk up the parent chain
      models_to_disconnect = parent_hierarchy(model)
      
      puts "Auto mode: Will disconnect records from #{models_to_disconnect.count} model(s) in hierarchy:"
      models_to_disconnect.each do |m|
        puts "  - #{m.model_name}"
      end
      puts

      # Disconnect in reverse order (children first, then parents)
      models_to_disconnect.reverse.each do |m|
        disconnect_model(m)
      end

      puts "\nAuto-disconnect complete for hierarchy."
    end

    def parent_hierarchy(model)
      # Walk up the parent chain similar to ReparentModelAction#model_tree
      # but in the opposite direction (parents, not children)
      models = [model]
      current_model = model

      # Walk up the parent chain
      while parent_attr = find_parent_attribute(current_model)
        parent_model = current_model.parent_model
        break unless parent_model
        break if parent_model.model_name.to_s == 'project' # Stop at project level
        
        models << parent_model
        current_model = parent_model
      end

      models
    end
  end
end
