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
      needs_hash = Metis::DataBlock.where(md5_hash: Metis::DataBlock::TEMP_MATCH).order(:updated_at).all[0..10]
      puts "Found #{needs_hash.count} data blocks to be checksummed."
      needs_hash.each(&:compute_hash!)

      needs_archive = Metis::DataBlock.exclude(md5_hash: Metis::DataBlock::TEMP_MATCH).where(archive_id: nil).order(:updated_at).all[0..10]
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

  class Assimilate < Etna::Command
    usage '<project> <bucket> <path> <files> Add the indicated (non-tracked) files to Metis at the given path'

    def execute project_name, bucket_name, folder_path, *files
      folder_path = folder_path.sub(%r!^/!,'')

      bucket = Metis::Bucket.where(project_name: project_name, name: bucket_name).first

      if folder_path.empty?
        files.each do |file_path|
          Metis::Assimilation.new(file_path, bucket).execute
        end
      else
        metis_folder = Metis::Folder.from_path(bucket, folder_path).last
        unless metis_folder
          puts "No such folder #{folder_path}"
          exit
        end
        files.each do |file_path|
          Metis::Assimilation.assimilate_folder(metis_folder, file_path)
        end
      end
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end
end
