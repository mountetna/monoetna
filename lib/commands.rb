class Metis
  class Help < Etna::Command
    usage 'List this help'

    def execute
      puts 'Commands:'
      Metis.instance.commands.each do |name,cmd|
        puts cmd.usage
      end
    end

    def setup(config)
      Metis.instance.configure(config)
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected magma instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      Metis.instance.configure(config)
      Metis.instance.load_models
    end
  end

  class Schema < Etna::Command
    usage 'Show the current database schema.'

    def execute
      Metis.instance.db.tap do |db|
        db.extension(:schema_dumper)
        puts db.dump_schema_migration
      end
    end

    def setup(config)
      super
      Metis.instance.setup_db
    end
  end

  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'

    def execute(version=nil)
      Sequel.extension(:migration)
      db = Metis.instance.db

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
      Metis.instance.setup_db
    end
  end

  class Archive < Etna::Command
    usage 'Checksum and archive files.'

    def execute
      needs_hash = Metis::File.where(file_hash:nil).order(:updated_at).all[0..10]
      puts "Found #{needs_hash.count} files to be checksummed."
      needs_hash.each(&:compute_hash!)

      needs_archive = Metis::File.exclude(file_hash: nil).where(archive_id: nil).order(:updated_at).all[0..10]
      puts "Found #{needs_archive.count} files to be archived."
      needs_archive.each do |file|
        begin
          file.backup!
        rescue ArgumentError => e
          puts "Could not archive #{file.file_name}"
          next
        end
      end
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end
end
