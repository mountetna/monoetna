class Metis
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

  class MeasureFileCounts < Etna::Command
    def setup(config)
      super
      Metis.instance.load_models
      Metis.instance.setup_db
    end

    # TODO: Also collect file size, which would require us to add that to the data_blocks table rather than
    # reading it from file.
    def execute
      Metis::File.select_group(:bucket_id).select_more{count(:id)}.each do |row|
        count = row[:count]
        bucket_id = row[:bucket_id]
        bucket = Metis::Bucket.where(id: bucket_id).first
        next if bucket.nil?

        project_name = bucket.project_name
        bucket_name = bucket.name

        tags = {project_name: project_name, bucket_name: bucket_name}
        Yabeda.metis.file_count.set(tags, count)
      end
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
    string_flags << '--version'

    def execute(version: nil)
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

  class ChecksumFiles < Etna::Command
    usage 'Checksum files.'

    def execute
      needs_hash = Metis::DataBlock.where(md5_hash: Metis::DataBlock::TEMP_MATCH, removed: false).order(:updated_at).all[0..10]
      puts "Found #{needs_hash.count} data blocks to be checksummed."
      needs_hash.each(&:compute_hash!)
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end

  class Archive < Etna::Command
    usage 'Archive files.'

    def execute
      needs_archive = Metis::DataBlock.exclude(md5_hash: Metis::DataBlock::TEMP_MATCH).where(archive_id: nil, removed: false).order(:updated_at).all[0..10]
      puts "Found #{needs_archive.count} files to be archived."
      needs_archive.each do |data_block|
        begin
          data_block.backup!
        rescue ArgumentError => e
          puts "Could not archive #{data_block.description}"
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

  class CreateDb < Etna::Command
    usage '# create the initial database per config.yml'

    def execute
      if @no_db
        create_db if @no_db

        puts "Database is setup. Please run `bin/metis migrate #{@project_name}`."
      else
        puts "Database already exists."
      end
    end

    def create_db
      # Create the database only

      puts "Creating database #{@db_config[:database]}"
      %x{ PGPASSWORD=#{@db_config[:password]} createdb -w -U #{@db_config[:user]} #{@db_config[:database]} }

      Metis.instance.setup_db
    end

    def setup(config)
      super
      @db_config = Metis.instance.config(:db)
      begin
        Metis.instance.setup_db
      rescue Sequel::DatabaseConnectionError
        @no_db = true
      end
    end
  end

  class RemoveOrphanDataBlocks < Etna::Command
    usage '# remove unused (orphaned) data blocks'

    def execute
      zero_hash = 'd41d8cd98f00b204e9800998ecf8427e'

      used_data_block_ids = Metis::File.all().map { |file| file.data_block.id }.uniq
      orphaned_data_blocks = Metis::DataBlock.exclude(id: used_data_block_ids).exclude(removed: true).exclude(md5_hash: zero_hash).all
      Metis.instance.logger.info("Found #{orphaned_data_blocks.count} orphaned data blocks to be removed.") if orphaned_data_blocks.count > 0
      orphaned_data_blocks.each do |orphaned_data_block|
        begin
          md5_hash_to_remove = orphaned_data_block.md5_hash
          orphaned_data_block.remove!
          Metis.instance.logger.info("Removed data_block with hash #{md5_hash_to_remove}")
        rescue Error => e
          Metis.instance.logger.log_error(e)
        end
      end
    end

    def setup(config)
      super
      Metis.instance.setup_logger
      Metis.instance.load_models
    end
  end
end
