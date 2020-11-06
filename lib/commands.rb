require 'date'
require 'logger'
require 'rollbar'
require 'tempfile'
require_relative 'helpers'
require 'yaml'

# /bin/etna will confirm execution before running any command that includes this module.
module RequireConfirmation
end

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
      raise "You have multiple environments configured, please specify your environment by adding --environment #{@config.keys.join("|")}"
    else
      raise "You do not have a successfully configured environment, please run #{program_name} config set https://polyphemus.ucsf.edu"
    end
  end

  def set_environment(env)
    @environment = env.nil? ? nil : env.to_sym
  end

  class Config
    include Etna::CommandExecutor

    class Show < Etna::Command

      boolean_flags << '--all'

      def execute(all: false)
        if all
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
        polyphemus_client ||= Etna::Clients::Polyphemus.new(
            host: host,
            token: token,
            ignore_ssl: ignore_ssl)
        workflow = Etna::Clients::Polyphemus::SetConfigurationWorkflow.new(
            polyphemus_client: polyphemus_client,
            config_file: EtnaApp.config_file_path)
        config = workflow.update_configuration_file(ignore_ssl: ignore_ssl)
        logger.info("Updated #{config.environment} configuration from #{host}.")
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

      class Add < Etna::Command
        include WithEtnaClients
        include WithLogger
        include StrongConfirmation

        boolean_flags << '--execute'
        string_flags << '--file'
        string_flags << '--target-model'

        def execute(project_name, execute: false, target_model: 'project', file: "#{project_name}_models_#{target_model}_tree.csv")
          reset

          unless File.exists?(file)
            puts "File #{file} is being prepared from the #{project_name} project."
            puts "Copying models descending from #{target_model}..."
            prepare_template(file, project_name, target_model)
            puts
            puts "Done!  You can start editing the file #{file} now, and I will report validation errors here."
          end

          load_models_from_csv(file)

          while true
            if @models && @errors.empty?
              puts "File #{file} is well formatted and contains #{@models.model_keys.length} models to synchronize to #{environment} #{project_name}."

              if execute
                puts "Would you like to execute?"
                if confirm
                  workflow.synchronize_to_server(@models, project_name, target_model) do |update_action|
                    puts "Executing #{update_action.to_h}..."
                  end
                  puts "Success!"
                end

                return
              else
                puts "To commit, run \033[1;31m#{program_name} #{project_name} --file #{file} --target-model #{target_model} --execute\033[0m"
              end
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
          @models = nil
          @last_load = Time.at(0)
        end

        def load_models_from_csv(file)
          reset

          @last_load = File.stat(file).mtime
          @models = File.open(file, 'r') do |f|
            workflow.prepare_models_from_csv(f) do |err|
              @errors << err
            end
          end

          return if @errors.empty?

          puts "Input file #{file} is invalid:"
          @errors.each do |err|
            puts  "  * " + err.gsub("\n", "\n\t")
          end
        end

        def prepare_template(file, project_name, target_model)
          tf = Tempfile.new
          begin
            File.open(tf.path, 'wb') { |f| workflow.write_models_templats_csv(f, project_name, target_model) }
            FileUtils.cp(tf.path, file)
          ensure
            tf.close!
          end
        end
      end

      class Attributes
        include Etna::CommandExecutor

        class ValidateActions < Etna::Command
          include WithEtnaClients
          include WithLogger

          def execute(project_name, filepath)
            # Use the workflow to validate instead of the AttributeActionsValidator directly,
            #   because we also need to check the actions against the project, i.e.
            #   do the right models exist, do the attributes already exist, etc.
            attribute_actions_workflow = Etna::Clients::Magma::AttributeActionsFromJsonWorkflow.new(
                magma_client: magma_client,
                project_name: project_name,
                filepath: filepath)
            # If the workflow initialized, then no errors!
            puts "Attribute Actions JSON is well-formatted and is valid for project #{project_name}."
          end
        end

        class ExecuteActions < Etna::Command
          include WithEtnaClientsByEnvironment
          include WithLogger
          include RequireConfirmation

          def execute(project_name, filepath)
            attribute_actions_workflow = Etna::Clients::Magma::AttributeActionsFromJsonWorkflow.new(
                magma_client: magma_client,
                project_name: project_name,
                filepath: filepath)
            attribute_actions_workflow.run!
          end
        end

        class UpdateFromCsv < Etna::Command
          include WithEtnaClients
          include WithLogger
          include RequireConfirmation

          def magma_crud
            @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
                magma_client: magma_client,
                project_name: @project_name)
          end

          def execute(project_name, model_name, filepath)
            @project_name = project_name

            update_attributes_workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflowSingleModel.new(
                magma_crud: magma_crud,
                project_name: project_name,
                model_name: model_name,
                filepath: filepath)
            update_attributes_workflow.update_attributes
          end
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
