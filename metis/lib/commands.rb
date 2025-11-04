class Metis
  class Console < Etna::Command
    usage "Open a console with a connected magma instance."

    def execute
      require "irb"
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
      Metis::File.select_group(:bucket_id).select_more { count(:id) }.each do |row|
        count = row[:count]
        bucket_id = row[:bucket_id]
        bucket = Metis::Bucket.where(id: bucket_id).first
        next if bucket.nil?

        project_name = bucket.project_name
        bucket_name = bucket.name

        tags = { project_name: project_name, bucket_name: bucket_name }
        Yabeda.metis.file_count.set(tags, count)
      end
    end
  end

  class ManualUpload < Etna::Command
    string_flags << "--project_name"
    string_flags << "--bucket_name"
    string_flags << "--path"

    def execute(file, project_name:, bucket_name:, path:)
      bucket = Metis::Bucket.find(project_name: project_name, name: bucket_name)
      db = Metis::DataBlock.create_from(::File.basename(file), file, true)
      folder_id = Metis::Folder.mkdir_p(bucket, path, project_name, 'etna.agent@ucsf.edu|system').map(&:id).last

      Metis::File.create(
        project_name: project_name,
        file_name: ::File.basename(file),
        folder_id: folder_id,
        bucket: bucket,
        author: 'etna.agent@ucsf.edu|system',
        data_block: db
      )

      puts "Done!"
    end

    def setup(config)
      super
      Metis.instance.setup_db
      Metis.instance.load_models
    end
  end

  class Schema < Etna::Command
    usage "Show the current database schema."

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
    usage "Run migrations for the current environment."
    string_flags << "--version"

    def execute(version: nil)
      Sequel.extension(:migration)
      db = Metis.instance.db

      if version
        puts "Migrating to version #{version}"
        Sequel::Migrator.run(db, "db/migrations", target: version.to_i)
      else
        puts "Migrating to latest"
        Sequel::Migrator.run(db, "db/migrations")
      end
    end

    def setup(config)
      super
      Metis.instance.setup_db
    end
  end

  class ChecksumFiles < Etna::Command
    usage "Checksum files."

    def execute
      while true
        needs_hash = Metis::DataBlock.where(md5_hash: Metis::DataBlock::TEMP_MATCH, removed: false).order(:updated_at).all[0..10]
        count = needs_hash.count
        puts "Found #{count} data blocks to be checksummed."
        needs_hash.each(&:compute_hash!)
        break if count == 0
      end
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end

  class Archive < Etna::Command
    usage "Archive files."

    def execute
      needs_archive = Metis::DataBlock.exclude(md5_hash: Metis::DataBlock::TEMP_MATCH).where(archive_id: nil, removed: false).order(:updated_at).all[0..10]
      puts "Found #{needs_archive.count} files to be archived."
      needs_archive.each do |data_block|
        begin
          puts "Archiving data block with hash: #{data_block.md5_hash}."
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

  class CreateDb < Etna::Command
    usage "# create the initial database per config.yml"

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

  class GenerateThumbnails < Etna::Command
    usage "Generate image thumbnail files."

    def execute
      needs_thumbnail_check = Metis::DataBlock.exclude(md5_hash: Metis::DataBlock::TEMP_MATCH).where(has_thumbnail: nil, removed: false).order(:updated_at).all[0..10]

      puts "Found #{needs_thumbnail_check.count} data blocks."

      has_thumbnail, needs_thumbnail = needs_thumbnail_check.partition do |data_block|
        data_block.thumbnail_in_cache?
      end

      puts "#{needs_thumbnail.count} images require thumbnails to be generated."

      needs_thumbnail.each do |data_block|
        begin
          data_block.generate_thumbnail
        rescue Metis::ThumbnailError => e
          Metis.instance.logger.error("Could not generate thumbnail for #{data_block.md5_hash}")
          Metis.instance.logger.log_error(e)
          next
        end
      end
    end

    def setup(config)
      super
      Metis.instance.setup_logger
      Metis.instance.load_models
    end
  end

  class ResetThumbnailFlag < Etna::Command
    usage "Reset has_thumbnail to nil for all data blcoks, if we need to regenerate."

    def execute
      Metis.instance.db[:data_blocks].update(has_thumbnail: nil)
    end

    def setup(config)
      super
      Metis.instance.load_models
    end
  end

  class BackfillDataBlockLedger < Etna::Command
    usage "Usage: bin/metis backfill_data_block_ledger [<project_name>] <phase>
           # backfill the datablock ledger
           # -links phase requires project_name
           # -orphaned phase is system-wide and does not require project_name
           # Must specify exactly one phase: -links or -orphaned"

    def execute(*args)
      # Parse arguments manually
      project_name = nil
      phase = nil
      flags_seen = []

      args.each do |arg|
        case arg
        when '-links'
          flags_seen << '-links'
          phase = :links if phase.nil?
        when '-orphaned'
          flags_seen << '-orphaned'
          phase = :orphaned if phase.nil?
        else
          project_name = arg unless arg.start_with?('-')
        end
      end

      if flags_seen.length > 1
        puts "Error: Cannot specify multiple phase options"
        puts "Choose exactly one: -links or -orphaned"
        return
      end

      if phase.nil?
        puts "Error: Exactly one phase option is required"
        puts ""
        puts "Usage: bin/metis backfill_data_block_ledger [<project_name>] <phase>"
        puts ""
        puts "Phases:"
        puts "  -links    Backfill link events for existing files (requires project_name)"
        puts "  -orphaned Backfill orphaned datablocks (system-wide, no project_name needed)"
        puts ""
        puts "Examples:"
        puts "  bin/metis backfill_data_block_ledger athena -links"
        puts "  bin/metis backfill_data_block_ledger -orphaned"
        puts ""
        return
      end

      case phase
      when :links
        if project_name.nil? || project_name.empty?
          puts "Error: Project name is required for -links phase"
          puts ""
          puts "Usage: bin/metis backfill_data_block_ledger <project_name> -links"
          return
        end

        # Verify project exists
        project_files_count = Metis::File.where(project_name: project_name).count
        if project_files_count == 0
          puts "Warning: No files found for project '#{project_name}'"
          puts "Please verify the project name is correct."
          puts ""
          response = ask_user("Continue anyway? (y/n): ")
          return unless response.downcase == 'y'
        end

        puts "Starting datablock ledger backfill for project: #{project_name}"
        puts "=" * 70
        puts ""
        puts "Running: Backfill link events for existing files"
        puts ""
        backfill_links(project_name)
        puts ""
        puts "=" * 70
        puts "Backfill complete for project: #{project_name}"
        puts "=" * 70
      when :orphaned
        puts "Starting datablock ledger backfill (system-wide)"
        puts "=" * 70
        puts ""
        puts "Running: Backfill orphaned datablocks"
        puts ""
        backfill_orphaned
        puts ""
        puts "=" * 70
        puts "Backfill complete (system-wide)"
        puts "=" * 70
      end
    end

    private

    def ask_user(prompt)
      print prompt
      gets.chomp
    end

    def backfill_links(project_name)
      puts "PHASE 1: Backfilling 'link_file_to_datablock' events for existing files"
      puts ""
      
      total_files = Metis::File.where(project_name: project_name).count
      puts "Found #{total_files} files in project '#{project_name}' to process."
      
      backfilled = 0
      skipped = 0
      errors = 0
      batch_size = 1000
      
      Metis::File.where(project_name: project_name).each_slice(batch_size) do |file_batch|
        file_batch.each do |file|
          begin
            # Skip if already in ledger
            existing = Metis::DataBlockLedger.where(
              file_id: file.id,
              event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK
            ).first
            
            if existing
              skipped += 1
              puts "Skipping file #{file.id} (md5: #{file.data_block.md5_hash}) - already has link event"
              next
            end
            
            # Create a 'link_file_to_datablock' event for each existing file
            Metis::DataBlockLedger.create(
              project_name: file.project_name,
              md5_hash: file.data_block.md5_hash,
              file_path: file.file_path,
              file_id: file.id,
              data_block_id: file.data_block_id,
              event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK,
              triggered_by: 'system_backfill',
              size: file.data_block.size,
              bucket_name: file.bucket.name,
              created_at: file.created_at || DateTime.now
            )
            
            backfilled += 1
            
            if (backfilled + skipped) % 1000 == 0
              puts "Progress: #{backfilled + skipped}/#{total_files} files processed (#{((backfilled + skipped).to_f / total_files * 100).round(2)}%)"
            end
          rescue => e
            errors += 1
            puts "Error backfilling file #{file.id}: #{e.message}"
          end
        end
      end
      
      puts ""
      puts "Phase 1 Results:"
      puts "- Link file to datablock events created: #{backfilled}"
      puts "- Skipped (already in ledger): #{skipped}"
      puts "- Errors: #{errors}" if errors > 0
    end

    def backfill_orphaned
      puts "PHASE 2: Finding orphaned datablocks and creating 'unlink' events"
      puts ""
      puts "Querying database for orphaned datablocks system-wide..."
      puts "(This identifies datablocks that are not referenced by any files)"
      puts "When a file is deleted via the API, we remove Metis::File records, but Metis::DataBlock records remain."
      puts ""

      # Get all data blocks that are still associated with a file
      all_used_datablock_ids = Metis::File
        .all
        .map { |file| file.data_block.id }
        .uniq
      
      # Get all data blocks that exist and are not removed
      all_datablock_ids = Metis::DataBlock
        .select_map(:id)
        .uniq
      
      if all_datablock_ids.empty?
        puts "No datablocks found in the system."
        puts "No orphaned datablocks to process."
        return
      end
      
      puts "Found #{all_datablock_ids.length} datablocks in total."
      puts "Found #{all_used_datablock_ids.length} datablocks referenced by files."
      
      # Orphaned = datablocks that exist but are not referenced by ANY files (system-wide)
      orphaned_datablock_ids = all_datablock_ids - all_used_datablock_ids
      
      # Get the actual datablock objects
      orphaned_datablocks = Metis::DataBlock
        .where(id: orphaned_datablock_ids)
        .all
      
      puts "Found #{orphaned_datablocks.length} orphaned datablocks (not referenced by any files)."
      
      if orphaned_datablocks.empty?
        puts "No orphaned datablocks to process."
        return
      end
      
      puts ""
      unlink_created = 0
      skipped = 0
      errors = 0
      total_size = 0
      
      orphaned_datablocks.each do |datablock|
        begin
          # Check if unlink event already exists for this datablock
          existing_unlink = Metis::DataBlockLedger.where(
            data_block_id: datablock.id,
            event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK
          ).first
          
          if existing_unlink
            # Already has unlink event
            skipped += 1
            puts "Skipping datablock #{datablock.id} (md5: #{datablock.md5_hash}) - already has unlink event"
            next
          end
          
          # Create synthetic 'unlink_file_from_datablock' event for orphaned datablock
          Metis::DataBlockLedger.create(
            project_name: nil, # We don't know the project name for orphaned datablocks
            md5_hash: datablock.md5_hash,
            file_path: nil,  # No specific file (orphaned)
            file_id: nil,
            data_block_id: datablock.id,
            event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK,
            triggered_by: 'system_backfill',
            size: datablock.size,
            bucket_name: nil,
            created_at: datablock.updated_at || DateTime.now
          )
          
          unlink_created += 1
          total_size += datablock.size
          
          if unlink_created % 100 == 0
            puts "Progress: #{unlink_created}/#{orphaned_datablocks.length} unlink events created"
          end
        rescue => e
          errors += 1
          puts "Error processing orphaned datablock #{datablock.id}: #{e.message}"
        end
      end
      
      puts ""
      puts "Phase 2 Results:"
      puts "- Unlink file from datablock events created: #{unlink_created}"
      puts "- Skipped (already has unlink event): #{skipped}" if skipped > 0
      # Format size appropriately (bytes, KB, or MB)
      if total_size < 1024
        puts "- Total orphaned size: #{total_size} bytes"
      elsif total_size < 1024 * 1024
        puts "- Total orphaned size: #{(total_size / 1024.0).round(2)} KB"
      else
        puts "- Total orphaned size: #{(total_size / 1024.0 / 1024.0).round(2)} MB"
      end
      puts "- Errors: #{errors}" if errors > 0
      
      if unlink_created > 0
        puts ""
        puts "These orphaned datablocks can now be vacuumed using the vacuum_datablocks API endpoint."
      end
    end

    def setup(config)
      super
      Metis.instance.setup_logger
      Metis.instance.load_models
    end
  end
end

