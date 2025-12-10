class Metis
  class DataBlockLedger < Sequel::Model(:data_block_ledger)
    plugin :timestamps, create: :created_at, update: false
    plugin :validation_helpers

    CREATE_DATABLOCK = 'create_datablock'
    RESOLVE_DATABLOCK = 'resolve_datablock'
    REUSE_DATABLOCK = 'reuse_datablock'
    LINK_FILE_TO_DATABLOCK = 'link_file_to_datablock'
    UNLINK_FILE_FROM_DATABLOCK = 'unlink_file_from_datablock' # Orphaned datablocks
    REMOVE_DATABLOCK = 'remove_datablock'

    EVENT_TYPES = [CREATE_DATABLOCK, RESOLVE_DATABLOCK, REUSE_DATABLOCK, LINK_FILE_TO_DATABLOCK, UNLINK_FILE_FROM_DATABLOCK, REMOVE_DATABLOCK].freeze

    SYSTEM_BACKFILL = 'system_backfill'
    CHECKSUM_COMMAND= 'command_checksum'

    def validate
      super
      # Allow project_name to be nil for:
      # 1. system_backfill orphaned datablock unlink events (orphaned datablocks don't have a known project)
      # 2. REMOVE_DATABLOCK events for backfilled vacuum (backfilled datablocks don't have a known project)
      required_fields = [:md5_hash, :data_block_id, :event_type, :created_at]
      unless (triggered_by == SYSTEM_BACKFILL && event_type == UNLINK_FILE_FROM_DATABLOCK) ||
             event_type == REMOVE_DATABLOCK
        required_fields << :project_name
      end
      validates_presence required_fields
      validates_includes EVENT_TYPES, :event_type
      validates_format /^(temp-[a-f0-9]{32}|[a-f0-9]{32})$/i, :md5_hash
    end

    def self.log_create(file, datablock, user)
      create(
        project_name: file.project_name,
        md5_hash: datablock.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: datablock.id,
        event_type: CREATE_DATABLOCK,
        triggered_by: user_identifier(user),
        size: datablock.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log create event: #{e.message}")
    end

  def self.log_link(file, datablock, user)
      create(
        project_name: file.project_name,
        md5_hash: datablock.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: datablock.id,
        event_type: LINK_FILE_TO_DATABLOCK,
        triggered_by: user_identifier(user),
        size: datablock.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log link event: #{e.message}")
    end

    def self.log_resolve(file, datablock, event)
      create(
        project_name: file.project_name,
        md5_hash: datablock.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: datablock.id,
        event_type: RESOLVE_DATABLOCK,
        triggered_by: event,
        size: datablock.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log resolve event: #{e.message}")
    end

    def self.log_deduplicate(file, datablock, user)
      create(
        project_name: file.project_name,
        md5_hash: datablock.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: datablock.id,
        event_type: REUSE_DATABLOCK,
        triggered_by: user_identifier(user),
        size: datablock.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log deduplicate event: #{e.message}")
    end

    def self.log_unlink(file, datablock, user)
      create(
        project_name: file.project_name,
        md5_hash: datablock.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: datablock.id,
        event_type: UNLINK_FILE_FROM_DATABLOCK,
        triggered_by: user_identifier(user),
        size: datablock.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log unlink event: #{e.message}")
    end

    def self.log_vacuum(datablock, project_name, user)
      create(
        project_name: project_name,
        md5_hash: datablock.md5_hash,
        file_path: nil,
        file_id: nil,
        data_block_id: datablock.id,
        event_type: REMOVE_DATABLOCK,
        triggered_by: user_identifier(user),
        size: datablock.size,
        bucket_name: nil,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log vacuum event: #{e.message}")
    end

    def self.user_identifier(user)
      return 'system' unless user
      
      # If user is a User object, return email
      return user.email if user.respond_to?(:email)
      
      # If user is a string in "email|name" format (from display_name), extract email
      if user.is_a?(String) && user.include?('|')
        return user.split('|').first
      end
      
      # Otherwise return as string
      user.to_s
    end

    def self.find_orphaned_datablocks(project_name, include_projects: [])
      if project_name.nil?
        raise Etna::Error, "Project name is required to find orphaned datablocks"
      end

      # Get all datablock IDs that this project has linked
      linked_datablock_ids = where(project_name: project_name)
        .where(event_type: [LINK_FILE_TO_DATABLOCK, REUSE_DATABLOCK])
        .select_map(:data_block_id)
        .uniq

      # Get datablock IDs that have been unlinked (no longer in use)
      unlinked_datablock_ids = where(
        project_name: project_name,
        event_type: UNLINK_FILE_FROM_DATABLOCK
      ).select_map(:data_block_id)
        .uniq

      # Get datablock IDs currently in use by files in the current project (safety check)
      used_datablock_ids = Metis::File
        .where(project_name: project_name)
        .select_map(:data_block_id)
        .uniq

      # Get datablock IDs already vacuumed
      vacuumed_datablock_ids = where(
        project_name: project_name,
        event_type: REMOVE_DATABLOCK
      ).select_map(:data_block_id)
        .uniq

      # Orphaned = (linked by this project) and (unlinked) and not (currently in use) and not (already vacuumed)
      orphaned_ids = (linked_datablock_ids & unlinked_datablock_ids) - used_datablock_ids - vacuumed_datablock_ids

      # If include_projects is specified, filter out datablocks used by projects not in the allowed list
      unless include_projects.empty?
        # Allowed projects = current project + included projects
        allowed_projects = [project_name] + include_projects
        
        # Get datablock IDs used by projects NOT in the allowed list
        blocking_datablock_ids = Metis::File
          .exclude(project_name: allowed_projects)
          .select_map(:data_block_id)
          .uniq
        
        # Remove datablocks that are used by non-included projects
        orphaned_ids = orphaned_ids - blocking_datablock_ids
      else
        # Default behavior: check globally to prevent vacuuming datablocks used by ANY other project
        used_by_other_projects = Metis::File
          .exclude(project_name: project_name)
          .select_map(:data_block_id)
          .uniq
        
        orphaned_ids = orphaned_ids - used_by_other_projects
      end

      Metis::DataBlock
        .where(id: orphaned_ids)
        .exclude(removed: true)
        .all
    end

    def self.find_orphaned_datablocks_backfilled
      # Get datablock IDs that have unlink events with SYSTEM_BACKFILL (backfilled records)
      backfilled_datablock_ids = where(
        triggered_by: SYSTEM_BACKFILL,
        event_type: UNLINK_FILE_FROM_DATABLOCK
      ).select_map(:data_block_id).uniq
      
      # Exclude datablock IDs that have already been vacuumed (have REMOVE_DATABLOCK events)
      vacuumed_datablock_ids = where(
        event_type: REMOVE_DATABLOCK
      ).select_map(:data_block_id).uniq
      
      orphaned_ids = backfilled_datablock_ids - vacuumed_datablock_ids
      
      Metis::DataBlock
        .where(id: orphaned_ids)
        .exclude(removed: true) # a bit redundant - remove will be set to true by the vacuum command (but this is also called in compute hash if file doesnt have data)
        .all
    end

    def self.calculate_event_counts(project_name = nil, include_projects: [])
      # Count all event types
      # If project_name is provided, count events for that project (and optionally include_projects)
      # If project_name is nil (backfilled mode), count only SYSTEM_BACKFILL events
      event_counts = if project_name
        all_projects = [project_name] + include_projects
        where(project_name: all_projects)
          .select_map(:event_type)
          .tally
      else
        # Backfilled mode: only count SYSTEM_BACKFILL events
        where(triggered_by: SYSTEM_BACKFILL)
          .select_map(:event_type)
          .tally
      end
      
      # Initialize all event types to 0 if they don't exist
      EVENT_TYPES.each do |event_type|
        event_counts[event_type] ||= 0
      end
      
      event_counts
    end

    def self.calculate_project_breakdown(orphaned_datablocks, project_name, include_projects = [])
      # Calculate project breakdown: count datablocks per project
      # For each orphaned datablock, count it for each project that has linked/reused it
      all_projects = [project_name] + include_projects
      project_breakdown = {}
      all_projects.each do |proj|
        project_breakdown[proj.to_sym] = 0
      end
      
      orphaned_datablocks.each do |datablock|
        # For each project, check if this datablock has been linked/reused in that project
        # If it has, increment that project's counter
        all_projects.each do |proj|
          # Check if this datablock has any link or reuse events in this project
          has_events = where(project_name: proj, data_block_id: datablock.id)
            .where(event_type: [LINK_FILE_TO_DATABLOCK, REUSE_DATABLOCK])
            .count > 0
          
          # If the datablock has events in this project, count it for this project
          # This means if a datablock is in both athena and labors, both get incremented
          project_breakdown[proj.to_sym] += 1 if has_events
        end
      end
      
      project_breakdown
    end

    def self.build_vacuum_details(orphaned_datablocks, project_name = nil)
      orphaned_datablocks.map do |datablock|
        # Get all events for this datablock to find any available file information
        all_events = if project_name
          where(project_name: project_name, data_block_id: datablock.id).all
        else
          where(data_block_id: datablock.id).all
        end
        
        # Try to get project_name from any event that has it (for backfilled mode)
        event_project_name = all_events.find { |e| e.project_name }&.project_name
        
        # Collect all unique file paths and bucket names from all events
        files = all_events.map do |event|
          next if event.file_path.nil? && event.bucket_name.nil?
          {
            file_path: event.file_path,
            bucket_name: event.bucket_name
          }
        end.compact.uniq
        
        # If no file information available from ledger, note the datablock description
        # which may contain the original filename
        description = datablock.description if files.empty? && datablock.description
        
        {
          data_block_id: datablock.id,
          md5_hash: datablock.md5_hash,
          size: datablock.size,
          description: description,
          files: files
        }
      end
    end
  end
end

