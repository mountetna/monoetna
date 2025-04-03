require_relative 'commands'
require_relative 'storage'
require_relative 'orchestration'
require_relative 'dependency_manager'
require_relative 'ssh_connection_pool'

class Vulcan
  include Etna::Application
  attr_reader :db, :ssh_pool

  def setup_ssh()
    Vulcan.instance.logger.info("Setting up SSH connection pool...")
    c = config(:ssh)
    
    begin
      if c[:proxy]
        Vulcan.instance.logger.info("Proxy configuration detected - setting up SSH connection with jump host...")
        Vulcan.instance.logger.info("Connecting through proxy host: #{c[:proxy][:host]}")
        @ssh_pool ||= SSHConnectionPool.new(
          c[:host], 
          c[:username], 
          c[:password],
          proxy: c[:proxy]
        )
      else
        Vulcan.instance.logger.info("Setting up direct SSH connection to #{c[:host]}...")
        @ssh_pool ||= SSHConnectionPool.new(c[:host], c[:username], c[:password])
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

  def setup_yabeda
    Yabeda.configure do
      group :vulcan do
        histogram :job_runtime do
          comment "Time spent by each cell, including storage and docker execution."
          unit :seconds
          tags [:script_hash]
          buckets [1, 5, 15, 60, 150, 300]
        end

        gauge :storage_disk_usage do
          comment "Amount (in bytes) used by vulcan storage directories"
          tags [:dir]
        end
      end

      collect do
        output = `du #{Vulcan.instance.config(:data_folder)} --max-depth 1` rescue ""
        output.split("\n").each do |line|
          parts = line.split("\t")
          if parts.length > 1
            bytes, dir = parts
            vulcan.storage_disk_usage.set({dir: dir}, bytes.to_i)
          end
        end
      end
    end

    super
  end

  def dependency_manager
    @dependency_manager ||= Vulcan::DependencyManager.new
  end
end
