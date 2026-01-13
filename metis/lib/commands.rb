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
    usage "Usage: bin/metis backfill_data_block_ledger [--project_name <project_name>] [--links] [--orphaned]
           # backfill the datablock ledger
           # --links phase requires --project_name
           # --orphaned phase is system-wide and does not require --project_name
           # Must specify exactly one phase: --links or --orphaned"

    boolean_flags << "--links"
    boolean_flags << "--orphaned"
    string_flags << "--project_name"

    def execute(links: false, orphaned: false, project_name: nil)
      # Validate that exactly one phase is specified
      if links && orphaned
        puts "Error: Cannot specify multiple phase options"
        puts "Choose exactly one: --links or --orphaned"
        return
      end

      if !links && !orphaned
        puts "Error: Exactly one phase option is required"
        puts ""
        puts "Usage: bin/metis backfill_data_block_ledger [--project_name <project_name>] [--links] [--orphaned]"
        puts ""
        puts "Phases:"
        puts "  --links    Backfill link events for existing files (requires --project_name)"
        puts "  --orphaned Backfill orphaned datablocks (system-wide, no --project_name needed)"
        puts ""
        puts "Examples:"
        puts "  bin/metis backfill_data_block_ledger --project_name athena --links"
        puts "  bin/metis backfill_data_block_ledger --orphaned"
        puts ""
        return
      end

      if links
        if project_name.nil? || project_name.empty?
          puts "Error: Project name is required for --links phase"
          puts ""
          puts "Usage: bin/metis backfill_data_block_ledger --project_name <project_name> --links"
          return
        end

        # Verify project exists
        project_files_count = Metis::File.where(project_name: project_name).count
        if project_files_count == 0
          puts "Error: No files found for project '#{project_name}'"
          puts "Please verify the project name is correct."
          puts ""
          puts "To see all project names, run:"
          puts "  SELECT DISTINCT project_name FROM files ORDER BY project_name;"
          return
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
      elsif orphaned
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

    def setup(config)
      super
      Metis.instance.load_models
    end

    private

    def backfill_links(project_name)
      puts "PHASE 1: Backfilling 'link_file_to_datablock' events for existing files"
      puts ""
      
      total_files = Metis::File.where(project_name: project_name).count
      puts "Found #{total_files} files in project '#{project_name}' to process."
      
      backfilled = 0
      skipped = 0
      errors = 0
      batch_size = 1000
      
      Metis::File.where(project_name: project_name).eager(:data_block, :bucket).each_slice(batch_size) do |file_batch|
        # Bulk check for existing link events (avoids N+1 queries)
        file_ids = file_batch.map(&:id)
        existing_link_file_ids = Metis::DataBlockLedger
          .where(
            file_id: file_ids,
            event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK
          )
          .select_map(:file_id)
          .to_set
        
        file_batch.each do |file|
          begin
            # Skip if already in ledger
            if existing_link_file_ids.include?(file.id)
              skipped += 1
              next
            end
            
            # Create a 'link_file_to_datablock' event for each existing file
            # Don't use log_link here because when we run backfill commands we disable real time ledger logging
            Metis::DataBlockLedger.create(
              project_name: file.project_name,
              md5_hash: file.data_block.md5_hash,
              file_path: file.file_path,
              file_id: file.id,
              data_block_id: file.data_block_id,
              event_type: Metis::DataBlockLedger::LINK_FILE_TO_DATABLOCK,
              triggered_by: Metis::DataBlockLedger::SYSTEM_BACKFILL,
              size: file.data_block.size,
              bucket_name: file.bucket.name,
              created_at: file.created_at || DateTime.now
            )
            
            backfilled += 1
          rescue => e
            errors += 1
            puts "Error backfilling file #{file.id}: #{e.message}"
          end
        end
        
        # Progress report after each batch
        processed = backfilled + skipped
        if processed % 5000 == 0 || processed == total_files
          progress_pct = (processed.to_f / total_files * 100).round(2)
          puts "Progress: #{processed}/#{total_files} (#{progress_pct}%) | Created: #{backfilled} | Skipped: #{skipped}"
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
        .select_map(:data_block_id)
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
      
      # Bulk check for existing unlink events (avoids N+1 queries)
      puts "Checking for existing unlink events..."
      existing_unlink_ids = Metis::DataBlockLedger
        .where(
          data_block_id: orphaned_datablock_ids,
          event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK
        )
        .select_map(:data_block_id)
        .to_set
      puts "Found #{existing_unlink_ids.size} datablocks that already have unlink events."
      puts ""
      
      orphaned_datablocks.each do |datablock|
        begin
          # Check if unlink event already exists for this datablock
          if existing_unlink_ids.include?(datablock.id)
            # Already has unlink event
            skipped += 1
            next
          end
          
          # Create synthetic 'unlink_file_from_datablock' event for orphaned datablock
          # Don't use log_unlink here because when we run backfill commands we disable real time ledger logging
          Metis::DataBlockLedger.create(
            project_name: nil, # We don't know the project name for orphaned datablocks
            md5_hash: datablock.md5_hash,
            file_path: nil,  # No specific file (orphaned)
            file_id: nil,
            data_block_id: datablock.id,
            event_type: Metis::DataBlockLedger::UNLINK_FILE_FROM_DATABLOCK,
            triggered_by: Metis::DataBlockLedger::SYSTEM_BACKFILL,
            size: datablock.size,
            bucket_name: nil,
            created_at: datablock.updated_at || DateTime.now
          )
          
          unlink_created += 1
          total_size += datablock.size
          
          if (unlink_created + skipped) % 5000 == 0
            progress_pct = ((unlink_created + skipped).to_f / orphaned_datablocks.length * 100).round(2)
            puts "Progress: #{unlink_created + skipped}/#{orphaned_datablocks.length} (#{progress_pct}%) | Created: #{unlink_created} | Skipped: #{skipped}"
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
      # Format size appropriately
      if total_size < 1024
        puts "- Total orphaned size: #{total_size} bytes"
      elsif total_size < 1024 * 1024
        puts "- Total orphaned size: #{(total_size / 1024.0).round(2)} KB"
      elsif total_size < 1024 * 1024 * 1024
        puts "- Total orphaned size: #{(total_size / 1024.0 / 1024.0).round(2)} MB"
      else
        puts "- Total orphaned size: #{(total_size / 1024.0 / 1024.0 / 1024.0).round(2)} GB"
      end
      puts "- Errors: #{errors}" if errors > 0
      
      if unlink_created > 0
        puts ""
        puts "These orphaned datablocks can now be vacuumed using the vacuum_datablocks API endpoint."
      end
    end
  end

  class LedgerStats < Etna::Command
    usage "Usage: bin/metis ledger_stats [--project_name <project_name>] [--backfilled] [--include_projects project1,project2]
           # Query ledger statistics and vacuum candidates
           # --project_name: Query stats for a specific project (tracked mode)
           # --backfilled: Query stats for backfilled datablocks
           # --include_projects: Optional comma-separated list of projects to include in vacuum calculations"

    string_flags << "--project_name"
    boolean_flags << "--backfilled"
    string_flags << "--include_projects"

    def execute(project_name: nil, backfilled: false, include_projects: nil)
      if !project_name && !backfilled
        puts "Error: Must specify either --project_name or --backfilled"
        puts ""
        puts "Usage: bin/metis ledger_stats [--project_name <project_name>] [--backfilled]"
        puts ""
        puts "Examples:"
        puts "  bin/metis ledger_stats --project_name athena"
        puts "  bin/metis ledger_stats --project_name athena --include_projects labors,victory"
        puts "  bin/metis ledger_stats --backfilled"
        return
      end

      if project_name && backfilled
        puts "Error: Cannot specify both --project_name and --backfilled"
        return
      end

      puts "Querying ledger statistics..."
      puts "=" * 70
      puts ""

      # Parse include_projects
      include_projects_array = include_projects ? include_projects.split(',').map(&:strip) : []

      # Use the service directly
      service = Metis::LedgerStatsService.new(
        project_name: project_name,
        backfilled: backfilled,
        include_projects: include_projects_array
      )
      
      begin
        data = service.calculate_stats
      rescue => e
        puts "Error: #{e.message}"
        return
      end
      
      puts "PROJECT: #{data[:project_name]}"
      puts ""
      
      # Display event counts
      puts "EVENT COUNTS:"
      puts "-" * 70
      event_counts = data[:event_counts] || {}
      event_counts.each do |event_type, count|
        formatted_event = event_type.to_s.split('_').map(&:capitalize).join(' ')
        puts "  #{formatted_event.ljust(30)} #{count.to_s.rjust(10)}"
      end
      puts ""
      
      # Display vacuum stats
      vacuum = data[:vacuum] || {}
      puts "VACUUM CANDIDATES:"
      puts "-" * 70
      puts "  Datablocks that can be vacuumed: #{vacuum[:datablocks_can_vacuum]}"
      puts "  Space that can be cleared:       #{format_size(vacuum[:space_can_clear])}"
      
      if vacuum[:include_projects] && !vacuum[:include_projects].empty?
        puts "  Include projects:                #{vacuum[:include_projects].join(', ')}"
      end
      
      if vacuum[:project_breakdown]
        puts ""
        puts "  Project Breakdown:"
        vacuum[:project_breakdown].each do |proj, count|
          puts "    #{proj.to_s.ljust(25)} #{count}"
        end
      end
      
      puts ""
      
      if vacuum[:details] && !vacuum[:details].empty?
        puts "VACUUM DETAILS (first 10 datablocks):"
        puts "-" * 70
        vacuum[:details].first(10).each do |detail|
          puts "  MD5: #{detail[:md5_hash]}"
          puts "  Size: #{format_size(detail[:size])}"
          puts "  Datablock ID: #{detail[:data_block_id]}"
          if detail[:description]
            puts "  Description: #{detail[:description]}"
          end
          if detail[:files] && !detail[:files].empty?
            puts "  Files:"
            detail[:files].first(3).each do |file|
              puts "    - #{file[:bucket_name]}/#{file[:file_path]}" if file[:file_path]
            end
          end
          puts ""
        end
        
        if vacuum[:details].length > 10
          puts "  ... and #{vacuum[:details].length - 10} more datablocks"
          puts ""
        end
      end
      
      puts "=" * 70
      if vacuum[:datablocks_can_vacuum] > 0
        puts "To vacuum these datablocks, run:"
        if backfilled
          puts "  bin/metis vacuum_datablocks backfilled --commit"
        else
          cmd = "  bin/metis vacuum_datablocks #{project_name}"
          cmd += " --include_projects #{include_projects}" if include_projects
          cmd += " --commit"
          puts cmd
        end
      else
        puts "No datablocks available to vacuum."
      end
      puts "=" * 70
    end

    def setup(config)
      super
      Metis.instance.load_models
    end

    private

    def format_size(bytes)
      return "0 bytes" if bytes.nil? || bytes == 0
      
      if bytes < 1024
        "#{bytes} bytes"
      elsif bytes < 1024 * 1024
        "#{(bytes / 1024.0).round(2)} KB"
      elsif bytes < 1024 * 1024 * 1024
        "#{(bytes / 1024.0 / 1024.0).round(2)} MB"
      else
        "#{(bytes / 1024.0 / 1024.0 / 1024.0).round(2)} GB"
      end
    end
  end

  class VacuumDatablocks < Etna::Command
    usage "Usage: bin/metis vacuum_datablocks <project_name> [--commit] [--include_projects project1,project2]
           # Vacuum (delete) orphaned datablocks
           # <project_name>: Project name or 'backfilled' for backfilled datablocks
           # --commit: Actually delete datablocks (defaults to dry-run)
           # --include_projects: Optional comma-separated list of projects to include"

    boolean_flags << "--commit"
    string_flags << "--include_projects"

    def execute(project_name, commit: false, include_projects: nil)
      if project_name.nil? || project_name.empty?
        puts "Error: Project name is required"
        puts ""
        puts "Usage: bin/metis vacuum_datablocks <project_name> [--commit] [--include_projects project1,project2]"
        puts ""
        puts "Examples:"
        puts "  bin/metis vacuum_datablocks athena                    # Dry-run"
        puts "  bin/metis vacuum_datablocks athena --commit           # Actually delete"
        puts "  bin/metis vacuum_datablocks backfilled --commit       # Vacuum backfilled datablocks"
        puts "  bin/metis vacuum_datablocks athena --include_projects labors --commit"
        return
      end

      if commit
        puts "WARNING: Running in COMMIT mode - datablocks will be DELETED!"
        puts "=" * 70
        puts ""
        print "Type 'yes' to confirm deletion: "
        confirmation = STDIN.gets.chomp
        if confirmation.downcase != 'yes'
          puts "Aborted."
          return
        end
        puts ""
      else
        puts "Running in DRY-RUN mode (no deletions will occur)"
        puts "Add --commit flag to actually delete datablocks"
        puts "=" * 70
        puts ""
      end

      # Parse include_projects
      include_projects_array = include_projects ? include_projects.split(',').map(&:strip) : []

      puts "Vacuuming datablocks for project: #{project_name}"
      puts ""

      # Use the service directly
      service = Metis::VacuumService.new(
        project_name: project_name,
        commit: commit,
        include_projects: include_projects_array,
        user: nil  # Commands run as system
      )
      
      begin
        data = service.vacuum_datablocks
      rescue => e
        puts "Error: #{e.message}"
        return
      end
      
      # Display results
      summary = data[:summary] || {}
      vacuumed = data[:vacuumed] || []
      errors = data[:errors] || []
      
      puts "VACUUM RESULTS:"
      puts "-" * 70
      puts "  Mode:                #{data[:dry_run] ? 'DRY-RUN' : 'COMMIT'}"
      puts "  Project:             #{summary[:project_name]}"
      if summary[:include_projects] && !summary[:include_projects].empty?
        puts "  Include projects:    #{summary[:include_projects].join(', ')}"
      end
      puts "  Total vacuumed:      #{summary[:total_vacuumed]}"
      puts "  Space freed:         #{format_size(summary[:space_freed])}"
      puts "  Errors:              #{summary[:errors_count]}"
      puts ""
      
      if vacuumed.length > 0
        puts "VACUUMED DATABLOCKS (first 20):"
        puts "-" * 70
        vacuumed.first(20).each_with_index do |db, idx|
          puts "  #{idx + 1}. MD5: #{db[:md5_hash]}"
          puts "     Size: #{format_size(db[:size])}"
          puts "     Location: #{db[:location]}" if db[:location]
          puts ""
        end
        
        if vacuumed.length > 20
          puts "  ... and #{vacuumed.length - 20} more datablocks"
          puts ""
        end
      end
      
      if errors.length > 0
        puts "ERRORS:"
        puts "-" * 70
        errors.each do |error|
          puts "  - #{error}"
        end
        puts ""
      end
      
      puts "=" * 70
      
      if data[:dry_run]
        puts "This was a DRY-RUN. No datablocks were actually deleted."
        puts "To commit the vacuum, run with --commit flag:"
        cmd = "  bin/metis vacuum_datablocks #{project_name}"
        cmd += " --include_projects #{include_projects}" if include_projects
        cmd += " --commit"
        puts cmd
      else
        puts "✓ Vacuum complete! #{summary[:total_vacuumed]} datablocks deleted."
        puts "✓ #{format_size(summary[:space_freed])} of space freed."
      end
      
      puts "=" * 70
    end

    def setup(config)
      super
      Metis.instance.load_models
    end

    private

    def format_size(bytes)
      return "0 bytes" if bytes.nil? || bytes == 0
      
      if bytes < 1024
        "#{bytes} bytes"
      elsif bytes < 1024 * 1024
        "#{(bytes / 1024.0).round(2)} KB"
      elsif bytes < 1024 * 1024 * 1024
        "#{(bytes / 1024.0 / 1024.0).round(2)} MB"
      else
        "#{(bytes / 1024.0 / 1024.0 / 1024.0).round(2)} GB"
      end
    end
  end
end

