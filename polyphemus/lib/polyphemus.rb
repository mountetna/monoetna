require_relative 'commands'
require_relative 'etl_cursor'
require 'etna/clients/magma'
require 'etna/clients/metis'

class Polyphemus
  include Etna::Application
  attr_reader :db

  def setup_db
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.extension :pg_json
    @db.pool.connection_validation_timeout = -1
  end
end
