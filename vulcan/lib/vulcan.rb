require_relative 'commands'
require_relative 'ssh_connection_pool'

class Vulcan
  include Etna::Application
  attr_reader :db, :ssh_pool

  def setup_ssh()
    Vulcan.instance.logger.info("Setting up SSH connection pool...")
    c = config(:ssh)
    begin
      password = c.key?(:password) ? c[:password] : nil
      use_ssh_config = c.key?(:use_ssh_config) ? c[:use_ssh_config] : false
      settings = c.key?(:settings) ? c[:settings] : {}
      @ssh_pool ||= SSHConnectionPool.new(c[:host], c[:username], password, use_ssh_config, settings: settings)
      @ssh_pool.with_conn do |ssh|
        result = ssh.exec!("echo 'SSH connection test successful'")
        if result.nil? || !result.include?('successful')
          raise "SSH test connection failed..."
        end
      end
    rescue => e
      Vulcan.instance.logger.error("Failed to establish SSH connection: #{e.message}")
      Vulcan.instance.logger.error(e.backtrace.join("\n"))
      raise e
    end
    
    Vulcan.instance.logger.info("SSH connection pool setup complete.")
  end

  def setup_db(load_models = true)
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.extension :pg_json
    @db.pool.connection_validation_timeout = -1

    require_relative 'models' if load_models
  end

  def vulcan_checks
    @remote_manager = Vulcan::RemoteManager.new(Vulcan.instance.ssh_pool)
    if config(:conda_env)
      Vulcan.instance.logger.info("Vulcan conda env: #{config(:conda_env)}")
      Vulcan.instance.logger.info("Attempting to activate conda environment...")
      command = @remote_manager.build_command
        .add('conda', 'activate', config(:conda_env))
        .add('snakemake', '--version')
      Vulcan.instance.logger.info("Running command: #{command.to_s}")
      result = @remote_manager.invoke_ssh_command(command.to_s, retries: 3)
      if result[:exit_status] == 0
        Vulcan.instance.logger.info("Snakemake version: #{result[:stdout].strip}")
      else
        Vulcan.instance.logger.error("Failed to activate conda environment and invoke Snakemake binary")
        raise "Snakemake binary does not exist or failed to execute"
      end
    else
      Vulcan.instance.logger.error("Please specify a conda environment in the config...")
    end
  end
end
