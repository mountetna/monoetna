require 'date'
require 'logger'
require 'readline'

class Vulcan
  class Schema < Etna::Command
    usage 'Show the current database schema.'

    def execute
      Vulcan.instance.db.tap do |db|
        db.extension(:schema_dumper)
        puts db.dump_schema_migration
      end
    end

    def setup(config)
      super
      Vulcan.instance.setup_db
    end
  end

  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'
    string_flags << '--version'

    def execute(version: nil)
      Sequel.extension(:migration)
      db = Vulcan.instance.db

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
      Vulcan.instance.setup_db(false)
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Vulcan instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
      Vulcan.instance.setup_db
    end
  end

  class CreateDb < Etna::Command
    usage '# create the initial database per config.yml'

    def execute
      if @no_db
        create_db if @no_db

        puts "Database is setup. Please run `bin/timur migrate #{@project_name}`."
      else
        puts "Database already exists."
      end
    end

    def create_db
      # Create the database only

      # We can't set environment vars for the executing shell,
      #   but we can set env vars for subshells that we open
      # https://stackoverflow.com/a/2660833
      fork do
        puts "Creating database #{@db_config[:database]}"
        ENV['PGPASSWORD'] = @db_config[:password]
        ENV['PGUSER'] = @db_config[:user]
        ENV['PGDB'] = @db_config[:database]
        %x{ createdb -w -U "$PGUSER" "$PGDB" }

        Vulcan.instance.setup_db
        exit 0
      end
      Process.wait
    end

    def setup(config)
      super
      @db_config = Vulcan.instance.config(:db)
      begin
        Vulcan.instance.setup_db
      rescue Sequel::DatabaseConnectionError
        @no_db = true
      end
    end
  end
end
