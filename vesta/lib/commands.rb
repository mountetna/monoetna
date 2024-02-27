require 'date'
require 'logger'
require 'rollbar'
require 'sequel'
require 'tempfile'
require 'active_support/all'
# require_relative 'helpers'


class Vesta
  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'
    string_flags << '--version'

    def execute(version: nil)
      Sequel.extension(:migration)
      db = Vesta.instance.db

      if version
        puts "Migrating to version #{version}"
        Sequel::Migrator.run(db, 'db/migrations', target: version.to_i)
      else
        puts 'Migrating to latest'
        Sequel::Migrator.run(db, 'db/migrations')
      end
    end

    def setup(config)
      super
      Vesta.instance.setup_db
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Vesta instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
      Vesta.instance.setup_db
      Vesta.instance.setup_sequel
      Vesta.instance.setup_ssh
    end
  end
end
