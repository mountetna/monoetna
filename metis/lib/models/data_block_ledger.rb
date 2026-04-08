class Metis
  class DataBlockLedger < Sequel::Model(:data_block_ledger)
    plugin :timestamps, create: :created_at, update: false
    plugin :validation_helpers

    CREATE_DATABLOCK = 'create_datablock'
    RESOLVE_DATABLOCK = 'resolve_datablock'
    REUSE_DATABLOCK = 'reuse_datablock'
    RESTORE_DATABLOCK = 'restore_datablock'
    LINK_FILE_TO_DATABLOCK = 'link_file_to_datablock'
    UNLINK_FILE_FROM_DATABLOCK = 'unlink_file_from_datablock' # Orphaned datablocks
    REMOVE_DATABLOCK = 'remove_datablock'

    EVENT_TYPES = [CREATE_DATABLOCK, RESOLVE_DATABLOCK, REUSE_DATABLOCK, RESTORE_DATABLOCK, LINK_FILE_TO_DATABLOCK, UNLINK_FILE_FROM_DATABLOCK, REMOVE_DATABLOCK].freeze

    SYSTEM_BACKFILL = 'system_backfill'
    CHECKSUM_COMMAND = 'command_checksum'
    CANCEL_UPLOAD = 'cancel_upload'

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

    def self.log_event(event_type:, datablock:, triggered_by:, project_name: nil, file_path: nil, file_id: nil, bucket_name: nil, event_meta: nil)
      create(
        project_name: project_name,
        md5_hash: datablock.md5_hash,
        file_path: file_path,
        file_id: file_id,
        data_block_id: datablock.id,
        event_type: event_type,
        triggered_by: user_identifier(triggered_by),
        size: datablock.size,
        bucket_name: bucket_name,
        event_meta: event_meta,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log #{event_type} event: #{e.message}")
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

    def self.find_orphaned_datablocks_for_vacuum(project_name)
      if project_name.nil?
        raise Etna::Error, "Project name is required to find orphaned datablocks"
      end

      orphaned_by_project, used_by_other_projects = orphaned_datablock_ids_for_project(project_name)

      # Block if ANY other project has a live file — vacuum is only safe when no one points to the datablock
      Metis::DataBlock.exclude_removed_and_temp_blocks(orphaned_by_project - used_by_other_projects)
    end

    # Returns datablocks this project has orphaned but cannot yet vacuum because other
    # projects still have live files pointing to them. Also returns which projects are
    # blocking each datablock so they can be tracked down.
    def self.find_blocked_orphaned_datablocks(project_name)
      if project_name.nil?
        raise Etna::Error, "Project name is required"
      end

      orphaned_by_project, used_by_other_projects = orphaned_datablock_ids_for_project(project_name)

      # Blocked = orphaned by this project AND still live in at least one other project
      blocked_ids = orphaned_by_project & used_by_other_projects

      datablocks = Metis::DataBlock.exclude_removed_and_temp_blocks(blocked_ids)

      # Build a map of datablock_id -> [blocking project names]
      blocking_files = Metis::File
        .exclude(project_name: project_name)
        .where(data_block_id: blocked_ids)
        .select(:data_block_id, :project_name)
        .all

      blocked_by = blocking_files.each_with_object({}) do |file, h|
        (h[file.data_block_id] ||= []) << file.project_name unless h[file.data_block_id]&.include?(file.project_name)
      end

      { datablocks: datablocks, blocked_by: blocked_by }
    end

    def self.orphaned_datablock_ids_for_project(project_name)
      linked_datablock_ids = where(project_name: project_name)
        .where(event_type: [LINK_FILE_TO_DATABLOCK, REUSE_DATABLOCK, RESTORE_DATABLOCK])
        .select_map(:data_block_id)
        .uniq

      unlinked_datablock_ids = where(
        project_name: project_name,
        event_type: UNLINK_FILE_FROM_DATABLOCK
      ).select_map(:data_block_id)
        .uniq

      live_datablock_ids = Metis::File
        .where(project_name: project_name)
        .select_map(:data_block_id)
        .uniq

      vacuumed_datablock_ids = where(
        project_name: project_name,
        event_type: REMOVE_DATABLOCK
      ).select_map(:data_block_id)
        .uniq

      restored_datablock_ids = where(
        project_name: project_name,
        event_type: RESTORE_DATABLOCK
      ).select_map(:data_block_id)
        .uniq

      # A block that was vacuumed but later restored is eligible again
      already_removed_ids = vacuumed_datablock_ids - restored_datablock_ids

      orphaned_by_project = (linked_datablock_ids & unlinked_datablock_ids) - live_datablock_ids - already_removed_ids

      used_by_other_projects = Metis::File
        .exclude(project_name: project_name)
        .select_map(:data_block_id)
        .uniq

      [orphaned_by_project, used_by_other_projects]
    end

    def self.find_orphaned_datablocks_backfilled_for_stats
      candidate_ids = backfilled_orphan_candidate_ids
      Metis::DataBlock.exclude_removed_and_temp_blocks(candidate_ids)
    end

    def self.find_orphaned_datablocks_backfilled_for_vacuum
      candidate_ids = backfilled_orphan_candidate_ids
      live_datablock_ids = Metis::File.select_map(:data_block_id).uniq

      # Block if any project has a live file — vacuum is only safe when no one points to the datablock
      Metis::DataBlock.exclude_removed_and_temp_blocks(candidate_ids - live_datablock_ids)
    end

    def self.backfilled_orphan_candidate_ids
      backfilled_datablock_ids = where(
        triggered_by: SYSTEM_BACKFILL,
        event_type: UNLINK_FILE_FROM_DATABLOCK
      ).select_map(:data_block_id).uniq

      vacuumed_datablock_ids = where(
        event_type: REMOVE_DATABLOCK
      ).select_map(:data_block_id).uniq

      restored_datablock_ids = where(
        event_type: RESTORE_DATABLOCK
      ).select_map(:data_block_id).uniq

      # A block that was vacuumed but later restored is eligible again
      already_removed_ids = vacuumed_datablock_ids - restored_datablock_ids

      backfilled_datablock_ids - already_removed_ids
    end

    def self.calculate_event_counts(project_name = nil)
      event_counts = if project_name
        where(project_name: project_name)
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

    def self.calculate_project_breakdown(orphaned_datablocks, project_name)
      project_breakdown = { project_name.to_sym => 0 }

      orphaned_datablocks.each do |datablock|
        has_events = where(project_name: project_name, data_block_id: datablock.id)
          .where(event_type: [LINK_FILE_TO_DATABLOCK, REUSE_DATABLOCK, RESTORE_DATABLOCK])
          .count > 0

        project_breakdown[project_name.to_sym] += 1 if has_events
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

