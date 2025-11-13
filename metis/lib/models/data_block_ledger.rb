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

    # Environment variable to enable/disable automatic logging
    # Set METIS_LEDGER_ENABLED=false to disable logging (useful when running backfill first)
    def self.ledger_enabled
      env_value = ENV['METIS_LEDGER_ENABLED']
      return false if env_value.nil?
      env_value.to_s.downcase == 'true'
    end

    def validate
      super
      # Allow project_name to be nil for:
      # 1. system_backfill orphaned datablock unlink events (orphaned datablocks don't have a known project)
      # 2. REMOVE_DATABLOCK events for legacy vacuum (legacy datablocks don't have a known project)
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
      return unless ledger_enabled
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
      return unless ledger_enabled
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
      return unless ledger_enabled
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
      return unless ledger_enabled
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
      return unless ledger_enabled
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
      return unless ledger_enabled
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
      user.respond_to?(:email) ? user.email : user.to_s
    end

    def self.find_orphaned_datablocks(project_name)
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

      # Get datablock IDs currently in use by files (safety check)
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

      Metis::DataBlock
        .where(id: orphaned_ids)
        .exclude(removed: true)
        .all
    end

    def self.find_orphaned_datablocks_legacy
      # Get datablock IDs that have unlink events with SYSTEM_BACKFILL (legacy backfilled records)
      backfilled_datablock_ids = where(
        triggered_by: SYSTEM_BACKFILL,
        event_type: UNLINK_FILE_FROM_DATABLOCK
      ).select_map(:data_block_id).uniq
      Metis::DataBlock
        .where(id: backfilled_datablock_ids)
        .all
    end
  end
end

