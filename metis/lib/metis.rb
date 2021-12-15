# metis.rb
require 'sequel'
require 'fileutils'
require_relative 'archiver'
require_relative 'assimilation'
require_relative 'thumbnail_cache'

# This class handles the http request and routing
class Metis
  include Etna::Application

  attr_reader :db

  def setup_db
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.pool.connection_validation_timeout = -1
  end

  def setup_logger
    # Override this to get daily logs instead of by size. Easier to
    #   debug Metis, since uploading in bulk won't overwrite
    #   errors.
    @logger = Etna::Logger.new(
      # The name of the log_file, required.
      config(:log_file),
      config(:log_age) || 'daily'  # Logger doesn't rotate out old daily logs, though, so we'll manage that via a system process
    )
    log_level = (config(:log_level) || 'warn').upcase.to_sym
    @logger.level = Logger.const_defined?(log_level) ? Logger.const_get(log_level) : Logger::WARN
  end

  def load_models
    setup_db

    require_relative 'models'
  end

  def project_path(project)
    ::File.join(config(:data_path), project)
  end

  def setup_yabeda
    Yabeda.configure do
      group :metis do
        gauge :file_count do
          comment "The number of files by project and bucket_name"
          tags [:project_name, :bucket_name]
        end
      end
    end

    super
  end

  def archiver
    @archiver ||= Metis::Archiver.new(config(:backup))
  end
end
