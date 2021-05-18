require_relative 'commands'
require_relative 'storage'
require_relative 'orchestration'

class Vulcan
  include Etna::Application
  attr_reader :db

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
          buckets [5, 15, 60, 150, 300]
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
end
