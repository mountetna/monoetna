# metis.rb
require 'sequel'
require 'fileutils'
require_relative 'archiver'

# This class handles the http request and routing
class Metis
  include Etna::Application

  attr_reader :db

  def setup_db
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.pool.connection_validation_timeout = -1
  end

  def load_models
    setup_db

    require_relative 'models'
  end

  def project_path(project)
    ::File.join(config(:data_path), project)
  end

  def archiver
    @archiver ||= Metis::Archiver.new(config(:backup))
  end
end
