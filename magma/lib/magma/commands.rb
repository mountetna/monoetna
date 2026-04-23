require 'date'
require 'logger'
require 'yaml'
require 'etna/command'
require_relative './actions/base_action'
require_relative './actions/set_model_template'

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

  class LinkTemplateModels < Etna::Command
    usage '[project_name] # Link DB-defined models to canonical template models'
    string_flags << '--template-project'
    string_flags << '--mapping-file'

    def execute(project_name, template_project: 'coprojects_template', mapping_file: default_mapping_file)
      template = Magma.instance.get_project(template_project)
      raise ArgumentError, "Template project '#{template_project}' does not exist" unless template

      mappings = read_template_mapping(mapping_file)
      link_project_models(project_name, template_project, template.models.values, mappings)
    end

    def setup(config)
      super
      Magma.instance.load_models
    end

    private

    def read_template_mapping(mapping_file)
      YAML.load_file(mapping_file)['mappings']
    end

    def link_project_models(project_name, template_project, template_models, mappings)
      project = Magma.instance.get_project(project_name)
      return unless project

      mappings.each do |mapping|
        template_model_name = mapping['template_model']
        aliases = mapping['potential_project_model_names'].map(&:downcase)

        template_model = template_models.find { |m| m.model_name.to_s == template_model_name }
        next unless template_model

        project.models.values.each do |model|
          next unless aliases.include?(model.model_name.to_s.downcase)
          next if model.template_model_name

          link_model(project_name, template_project, model, template_model)
        end
      end
    end

    def link_model(project_name, template_project, model, template_model)
      action = Magma::SetModelTemplateAction.new(project_name.to_s, {
        action_name: 'set_model_template',
        model_name: model.model_name.to_s,
        template_project_name: template_project,
        template_model_name: template_model.model_name.to_s,
      })

      if action.validate && action.perform
        puts "Linked #{project_name}.#{model.model_name} -> #{template_project}.#{template_model.model_name}"
      else
        puts "Failed to link #{project_name}.#{model.model_name}: #{action.errors.map { |e| e[:message] }.join(', ')}"
      end
    end

    def default_mapping_file
      File.expand_path('template_model_name_mappings.yml', __dir__)
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
end
