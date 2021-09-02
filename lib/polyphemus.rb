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

  def setup_ssh
    ssh_dir = ::File.expand_path(".ssh", "~")
    FileUtils.mkdir_p(ssh_dir)
    ::File.open(::File.join(ssh_dir, "known_hosts"), "w") do |known_hosts|
      Polyphemus.instance.config(:ingest).each do |ingest_config_type, ingest_config|
        next unless ingest_config.is_a?(Array)

        ingest_config.each do |conf|
          known_hosts.write(conf[:ssh_key] + "\n") unless conf[:ssh_key].nil?
        end
      end
    end
  end
end
