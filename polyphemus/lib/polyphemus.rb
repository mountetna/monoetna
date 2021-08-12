require_relative 'commands'
require_relative 'etl_cursor'
require_relative 'metis_file_etl'
require_relative 'magma_record_etl'
require_relative 'etls'
require_relative 'metis_file_for_magma_model_etl'
require 'etna/clients/magma'
require 'etna/clients/metis'
require 'etna/clients/janus'

class Polyphemus
  include Etna::Application
  attr_reader :db

  def setup_db
    return if @db
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.extension :pg_json
    @db.pool.connection_validation_timeout = -1
  end

  def setup_sequel
    Sequel::Model.plugin :timestamps, update_on_create: true
    Sequel::Model.require_valid_table = false
    Sequel.extension :inflector
    Sequel.extension :migration
    require_relative 'models/models'
  end
end
