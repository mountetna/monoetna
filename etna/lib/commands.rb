require 'date'
require 'logger'
require 'rollbar'
require 'tempfile'
require_relative 'helpers'
require 'yaml'

class EtnaApp
  def self.config_file_path
    File.join(Dir.home, 'etna.yml')
  end

  string_flags << '--environment'

  def dispatch_to_subcommand(cmd = 'help', *args, environment: nil, **kwds)
    set_environment(environment)
    super(cmd, *args, **kwds)
  end

  def environment
    if @environment
      @environment
    elsif @config && @config.is_a?(Hash) && @config.keys.length == 1
      @config.keys.last.to_sym
    elsif @config && @config.is_a?(Hash) && @config.keys.length > 1
      :many
    else
      :none
    end
  end

  def set_environment(env)
    @environment = env.nil? ? nil : env.to_sym
  end

  class Config
    include Etna::CommandExecutor

    class Show < Etna::Command
      def execute
        if EtnaApp.instance.environment == :many
          File.open(EtnaApp.config_file_path, 'r') { |f| puts f.read }
        else
          puts "Current environment: #{EtnaApp.instance.environment}"
          pp EtnaApp.instance.env_config
        end
      end
    end

    class Set < Etna::Command
      include WithEtnaClients
      include WithLogger

      boolean_flags << '--ignore-ssl'

      def execute(host, ignore_ssl: false)
        polyphemus_client = Etna::Clients::Polyphemus.new(
            host: host,
            token: token(ignore_environment: true),
            ignore_ssl: ignore_ssl)
        workflow = Etna::Clients::Polyphemus::SetConfigurationWorkflow.new(
            polyphemus_client: polyphemus_client,
            config_file: EtnaApp.config_file_path)
        config = workflow.update_configuration_file(ignore_ssl: ignore_ssl)
        logger&.info("Updated #{config.environment} configuration from #{host}.")
      end

      def setup(config)
        super
        EtnaApp.instance.setup_logger
      end
    end
  end

  class CreateTemplate
    include Etna::CommandExecutor

    # TODO: Refactor and replace with a command workflow similar to AddProjectModels
    class AttributeActions < Etna::Command
      def execute
        spec = Gem::Specification.find_by_name("etna")
        gem_root = spec.gem_dir
        FileUtils.cp(
            "#{gem_root}/lib/etna/templates/attribute_actions_template.json",
            'attribute_actions_template.json')
        puts "A sample attribute actions JSON template has been provided in the current directory as `attribute_actions_template.json`."
      end
    end
  end

  class Administrate
    include Etna::CommandExecutor

    class Token
      include Etna::CommandExecutor

      class Generate < Etna::Command
        include WithLogger

        boolean_flags << "--task"
        string_flags << "--project-name"
        string_flags << "--email"

        def execute(email:, task: false, project_name: nil)
          # the token is not required, but can be used if available
          # to generate a task token, so we pass it in here
          janus_client = Etna::Clients::Janus.new(
            token: ENV['TOKEN'],
            ignore_ssl: EtnaApp.instance.config(:ignore_ssl),
            **EtnaApp.instance.config(:janus, EtnaApp.instance.environment))

          generate_token_workflow = Etna::Clients::Janus::GenerateTokenWorkflow.new(
            janus_client: janus_client,
            token_type: task ? 'task' : 'login',
            email: email,
            project_name: project_name,
            private_key_file: EtnaApp.instance.config(:private_key, EtnaApp.instance.environment)
          )
          generate_token_workflow.generate!
        end
      end
    end

    class Project
      include Etna::CommandExecutor

      class Create < Etna::Command
        include WithEtnaClients
        include WithLogger

        def execute(project_name, project_name_full)
          create_args = {project_name: project_name, project_name_full: project_name_full}
          Etna::Clients::Magma::ProjectValidator.new(**create_args).validate!("#{program_name} args invalid!")

          create_project_workflow = Etna::Clients::Magma::CreateProjectWorkflow.new(
              magma_client: magma_client,
              janus_client: janus_client,
              **create_args
          )
          create_project_workflow.create!
        end
      end
    end

    class Models
      include Etna::CommandExecutor

      class CopyTemplate < Etna::Command
        include WithEtnaClients
        include WithLogger
        include StrongConfirmation

        string_flags << '--file'
        string_flags << '--target-model'

        def execute(project_name, target_model: 'project', file: "#{project_name}_models_#{target_model}_tree.csv")
          unless File.exists?(file)
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

        def prepare_template(file, project_name, target_model)
          tf = Tempfile.new

          begin
            File.open(tf.path, 'wb') { |f| workflow.write_models_template_csv(project_name, target_model, io: f) }
            FileUtils.cp(tf.path, file)
          ensure
            tf.close!
          end
        end
      end

      class ApplyTemplate < Etna::Command
        include WithEtnaClients
        include StrongConfirmation
        include WithLogger

        string_flags << '--file'
        string_flags << '--target-model'

        def execute(project_name, target_model: 'project', file: "#{project_name}_models_#{target_model}_tree.csv")
          reset

          unless File.exists?(file)
            puts "Could not find file #{file}"
            return
          end

          load_models_from_csv(file)

          while true
            if @changeset && @errors.empty?
              puts "File #{file} is well formatted.  Calculating expected changes..."
              sync_workflow = workflow.plan_synchronization(@changeset, project_name, target_model)
              models_and_action_types = sync_workflow.planned_actions.map { |a| Etna::Clients::Magma::ModelSynchronizationWorkflow.models_affected_by(a).map { |m| [m, a.action_name] }.flatten }
              models_and_action_types.sort!
              models_and_action_types = models_and_action_types.group_by(&:first)

              models_and_action_types.each do |model, actions|
                actions = actions.map { |a| a[1] }.sort
                actions = actions.group_by { |v| v }
                puts
                puts "#{model} changes:"
                actions.each do |type, actions|
                  puts " * #{type}: #{actions.length}"
                end
              end

              puts
              puts "Would you like to execute?"
              if confirm
                sync_workflow.update_block = Proc.new do |action|
                  puts "Executing #{action.action_name} on #{Etna::Clients::Magma::ModelSynchronizationWorkflow.models_affected_by(action)}..."
                end

                sync_workflow.execute_planned!
                File.unlink(file)
              end

              return
            end

            # Poll for updates
            puts "Watching for changes to #{file}..."
            while File.stat(file).mtime == @last_load
              sleep(1)
            end

            load_models_from_csv(file)
          end
        end

        def workflow
          @workflow ||= Etna::Clients::Magma::AddProjectModelsWorkflow.new(magma_client: magma_client)
        end

        def reset
          @errors = []
          @changeset = nil
          @last_load = Time.at(0)
        end

        def load_models_from_csv(file)
          reset

          @last_load = File.stat(file).mtime
          @changeset = File.open(file, 'r') do |f|
            workflow.prepare_changeset_from_csv(io: f) do |err|
              @errors << err
            end
          end

          return if @errors.empty?

          puts "Input file #{file} is invalid:"
          @errors.each do |err|
            puts "  * " + err.gsub("\n", "\n\t")
          end
        end
      end

      class Attributes
        include Etna::CommandExecutor

        class UpdateFromCsv < Etna::Command
          include WithEtnaClients
          include WithLogger

          boolean_flags << '--json-values'
          string_flags << '--hole-value'

          def magma_crud
            @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
                magma_client: magma_client,
                project_name: @project_name)
          end

          def execute(project_name, model_name, filepath, hole_value: '_', json_values: false)
            @project_name = project_name

            update_attributes_workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflowSingleModel.new(
                magma_crud: magma_crud,
                project_name: project_name,
                model_name: model_name,
                filepath: filepath,
                hole_value: hole_value,
                json_values: json_values)
            update_attributes_workflow.update_attributes
          end
        end

        class CreateFileLinkingCsv < Etna::Command
          include WithEtnaClients

          string_flags << '--file'
          string_flags << '--regex'
          string_flags << '--folder'
          boolean_flags << '--collection'

          def execute(project_name, bucket_name, attribute_name, extension, collection: false, regex: "**/*/(?<identifier>.+)\\.#{extension}$", folder: "", file: "#{project_name}_#{attribute_name}.csv")
            if folder.start_with?("/")
              folder = folder.slice(1..-1)
            end

            regex = Regexp.new(regex)

            workflow = Etna::Clients::Magma::SimpleFileLinkingWorkflow.new(
                metis_client: metis_client,
                project_name: project_name,
                bucket_name: bucket_name,
                folder: folder,
                extension: extension,
                attribute_name: attribute_name,
                regex: regex,
                file_collection: collection,
            )

            workflow.write_csv_io(filename: file)
          end
        end

        class LoadTableFromCsv < Etna::Command
          include WithEtnaClients

          boolean_flags << '--execute'

          def execute(project_name, model_name, file_path, execute: false)
            request = Etna::Clients::Magma::RetrievalRequest.new(project_name: project_name)
            request.model_name = model_name
            request.attribute_names = 'all'
            request.record_names = []
            model = magma_client.retrieve(request).models.model(model_name)
            model_parent_name = model.template.attributes.all.select do |attribute|
              attribute.attribute_type == Etna::Clients::Magma::AttributeType::PARENT
            end.first.name

            other_attribute_names = model.template.attributes.all.reject do |attribute|
              attribute.attribute_type == Etna::Clients::Magma::AttributeType::PARENT
            end.map do |attribute|
              attribute.name
            end

            # NOTE: This does not call ensure_parent currently because of MVIR1 consent--
            #   if the timepoint doesn't exist, the patient may be no study? (one example, at least)
            update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project_name)

            data = CSV.parse(File.read(file_path), headers: true)

            data.by_row.each do |row|
              revision = {}
              other_attribute_names.each do |attribute_name|
                revision[attribute_name] = row[attribute_name] unless row[attribute_name].nil?
              end
              update_request.append_table(model_parent_name, row[model_parent_name], model_name, revision)
            end

            puts update_request

            if execute
              magma_client.update_json(update_request)
            end
          end
        end
      end

      class LoadFromRedcap < Etna::Command
        include WithEtnaClients
        include WithLogger
        include StrongConfirmation

        boolean_flags << '--commit'
        string_flags << '--models'
        string_flags << '--mode'

        def execute(project_name, redcap_tokens, models: "all", mode: nil, commit: false)
          raise "Must provide at least one REDCap token (comma-separated)." unless redcap_tokens.split(',').length > 0

          puts "NOTE: This is a **preview** of what the data loading will look like. Use the --commit flag to load records into Magma." unless commit

          polyphemus_client.job(Etna::Clients::Polyphemus::RedcapJobRequest.new(
            model_names: "all" == models ? "all" : models.split(','),
            mode: mode,
            redcap_tokens: redcap_tokens.split(','),
            project_name: project_name,
            commit: commit
          )) do |response|
            response.read_body do |chunk|
              puts clean_sne_message(chunk)
            end
          end
        end

        def clean_sne_message(chunk)
          chunk.split("\n").reject do |c|
            c.start_with?("retry:") || c.start_with?("event:")
          end.map do |c|
            c.gsub("data:", "").strip
          end.reject do |c|
            c.empty?
          end
        end
      end

      class SetDateShiftRoot < Etna::Command
        include WithEtnaClients
        include WithLogger
        include StrongConfirmation

        boolean_flags << '--date-shift-root'
        string_flags << '--target-model'

        def execute(project_name, target_model: 'subject', date_shift_root: false)
          magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            actions:[Etna::Clients::Magma::SetDateShiftRootAction.new(
              model_name: target_model,
              date_shift_root: date_shift_root
            )],
            project_name: project_name,
          ))
        end
      end
    end
  end

  class Polyphemus < Etna::Command
    include WithEtnaClients
    string_flags << '--tag'
    boolean_flags << '--local'

    def execute(tag: nil, local: false)
      @tag = tag || default_tag
      @local = local
      create_config_yml

      unless @local
        puts `docker pull #{image}`
      end

      puts run_cmd
      exec(run_cmd)
    end

    def run_cmd
      parts = ['docker run --rm -it']

      unless Gem.win_platform?
        uid = Etc.getpwnam(ENV['USER']).uid
        gid = Etc.getgrnam('docker').gid
        parts << "-u #{uid}:#{gid}"
        parts << "-v /etc/passwd:/etc/passwd:ro"
        parts << "-v /etc/group:/etc/group:ro"
      end

      parts << "-v #{temp_config_file.path}:/app/config.yml:ro"
      parts << "-v #{Dir.pwd}:/app/workspace"
      parts << "-e HOME=/root"

      parts << image
      parts << "bash -c 'cd /app/workspace && exec bash --rcfile /root/.bashrc'"

      parts.join(' ')
    end

    def temp_config_file
      @temp_config_file ||= Tempfile.new
    end

    def create_config_yml
      host_config = [:metis, :magma, :janus, :timur].map do |host|
        [host, (EtnaApp.instance.config(host) || {}).update(token: token)]
      end.to_h

      config = {
          EtnaApp.instance.environment => {
                                              log_file: '/dev/stdout',
                                              log_level: 'info',
                                          }.update(host_config)
      }

      File.open(temp_config_file.path, 'w') { |f| YAML.dump(config, f) }
    end

    def default_tag
      docker = EtnaApp.instance.config(:docker)
      if docker.nil?
        nil
      else
        docker[:default_tag]
      end || 'master'
    end

    def image
      if @local
        "polyphemus:#{@tag}"
      else
        "etnaagent/polyphemus:#{@tag}"
      end
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Etna instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
    end
  end
end
