require_relative 'commands'
require_relative 'etl_cursor'
require_relative 'metis_file_etl'
require_relative 'magma_record_etl'
require_relative 'etls'
require 'etna/clients/magma'
require 'etna/clients/metis'
require 'etna/clients/janus'

class Polyphemus
  include Etna::Application
  attr_reader :db

  def setup_db
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.extension :pg_json
    @db.pool.connection_validation_timeout = -1
  end

  def setup_yabeda
    Yabeda.configure do
      group :polyphemus do
        gauge :last_command_completion do
          comment "Unix time of last time command was completed"
          tags [:command, :status]
        end
      end
    end

    super
  end

  def run_command(config, *args, &block)
    cmd, cmd_args, cmd_kwds = find_command(*args)

    begin
      super
      Yabeda.polyphemus.last_command_completion.set({ command: cmd.class.name, status: 'success' }, Time.now.to_i)
    rescue
      Yabeda.polyphemus.last_command_completion.set({ command: cmd.class.name, status: 'failed' }, Time.now.to_i)
      raise
    ensure
      write_job_metrics("#{cmd.class.name}.last_run")
    end
  end
end
