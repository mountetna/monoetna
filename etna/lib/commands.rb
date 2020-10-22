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
    elsif @config && @config.keys.length == 1
      @config.keys.last.to_sym
    else
      :production
    end
  end

  def set_environment(env)
    @environment = env.nil? ? nil : env.to_sym
  end

  class Config
    include Etna::CommandExecutor

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

  class Administrate
    include Etna::CommandExecutor

    class CreateProject < Etna::Command
      include WithEtnaClients
      include WithLogger

      def execute(filepath)
        create_project_workflow = Etna::Clients::Magma::CreateProjectFromJsonWorkflow.new(
            magma_client: magma_client,
            janus_client: janus_client,
            filepath: filepath
        )
        create_project_workflow.create!
      end
    end

    class ValidateProject < Etna::Command
      include WithLogger
      usage 'validate_project <filepath>'

      def execute(filepath)
        project = Etna::Clients::Magma::JsonProject.new(filepath: filepath)
        return puts "Project JSON is well-formatted!" if project.valid?

        puts "Project JSON has #{project.errors.length} errors:"
        project.errors.each do |error|
          puts "  * " + error.gsub("\n", "\n\t")
        end
      end
    end

    class ValidateModel < Etna::Command
      include WithEtnaClients
      include WithLogger

      def execute(project_name, model_name, filepath)
        Etna::Clients::Magma::AddModelFromJsonWorkflow.new(
            magma_client: magma_client,
            project_name: project_name,
            model_name: model_name,
            filepath: filepath)
        # If the workflow initialized, then no errors!
        puts "Model JSON for #{model_name} is well-formatted and is valid for project #{project_name}."
      end
    end

    class AddModel < Etna::Command
      include WithEtnaClients
      include WithLogger
      include RequireConfirmation

      def execute(project_name, model_name, filepath)
        add_model_workflow = Etna::Clients::Magma::AddModelFromJsonWorkflow.new(
            magma_client: magma_client,
            project_name: project_name,
            model_name: model_name,
            filepath: filepath)
        add_model_workflow.add!
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
    usage 'Open a console with the local etna gem environment'

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
