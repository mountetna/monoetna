require_relative 'commands'
require_relative 'etl_cursor'
require_relative 'metis_file_etl'
require_relative 'magma_record_etl'
require_relative 'etls'
require 'etna/clients/magma'
require 'etna/clients/metis'
require 'etna/clients/janus'
require 'pathname'

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

  def setup_yabeda
    Yabeda.configure do
      group :linker do
        gauge :linked_attributes do
          comment "The number of files linked for each attribute"
          tags [:project_name, :model_name, :attribute_type, :attribute_name]
        end

        gauge :identified_records do
          comment "The number of records identified by the linker"
          tags [:project_name, :model_name]
        end
      end
    end

    super
  end

  def setup_ssh
    ssh_dir = ::File.expand_path(".ssh", "~")
    FileUtils.mkdir_p(ssh_dir)

    known_hosts_file = ::File.join(ssh_dir, "known_hosts")

    ::File.open(known_hosts_file, "w") do |known_hosts|
      Polyphemus.instance.config(:ingest).each do |ingest_config_type, ingest_config|
        next unless ingest_config.is_a?(Array)

        ingest_config.each do |conf|
          known_hosts.write(conf[:ssh_key] + "\n") unless conf[:ssh_key].nil?
        end
      end
    end if Pathname.new(known_hosts_file).writable?
  end
end
