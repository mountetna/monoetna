require_relative 'commands'
require 'pathname'

class Gnomon
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
    #require_relative 'models/models'
  end
end
