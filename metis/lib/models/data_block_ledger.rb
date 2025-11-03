class Metis
  class DataBlockLedger < Sequel::Model(:data_block_ledger)
    plugin :timestamps, create: :created_at, update: false
    plugin :validation_helpers

    CREATE_DATABLOCK = 'create_datablock'
    REUSE_DATABLOCK = 'reuse_datablock'
    LINK_FILE_TO_DATABLOCK = 'link_file_to_datablock'
    UNLINK_FILE_FROM_DATABLOCK = 'unlink_file_from_datablock'
    REMOVE_DATABLOCK = 'remove_datablock'

    EVENT_TYPES = [CREATE_DATABLOCK, REUSE_DATABLOCK, LINK_FILE_TO_DATABLOCK, UNLINK_FILE_FROM_DATABLOCK, REMOVE_DATABLOCK].freeze

    SYSTEM_BACKFILL = 'system_backfill'

    def validate
      super
      # Allow project_name to be nil for system_backfill orphaned datablock unlink events
      # (orphaned datablocks don't have a known project since the file was deleted)
      required_fields = [:md5_hash, :data_block_id, :event_type, :created_at]
      unless triggered_by == SYSTEM_BACKFILL && event_type == UNLINK_FILE_FROM_DATABLOCK
        required_fields << :project_name
      end
      validates_presence required_fields
      validates_includes EVENT_TYPES, :event_type
      validates_format /^[a-f0-9]{32}$/i, :md5_hash
    end

    def self.log_create(datablock:, project_name:, user:)
      create(
        project_name: project_name,
        md5_hash: datablock.md5_hash,
        file_path: nil,
        file_id: nil,
        data_block_id: datablock.id,
        event_type: CREATE_DATABLOCK,
        triggered_by: user_identifier(user),
        size: datablock.size,
        bucket_name: nil,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log create event: #{e.message}")
    end

    def self.log_deduplicate(old_datablock:, new_datablock:, file:, user:)
      create(
        project_name: file.project_name,
        md5_hash: new_datablock.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: new_datablock.id,
        event_type: REUSE_DATABLOCK,
        triggered_by: user_identifier(user),
        size: new_datablock.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log deduplicate event: #{e.message}")
    end

    def self.log_link(file:, user:)
      create(
        project_name: file.project_name,
        md5_hash: file.data_block.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: file.data_block_id,
        event_type: LINK_FILE_TO_DATABLOCK,
        triggered_by: user_identifier(user),
        size: file.data_block.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log link event: #{e.message}")
    end

    def self.log_unlink(file:, user:)
      create(
        project_name: file.project_name,
        md5_hash: file.data_block.md5_hash,
        file_path: file.file_path,
        file_id: file.id,
        data_block_id: file.data_block_id,
        event_type: UNLINK_FILE_FROM_DATABLOCK,
        triggered_by: user_identifier(user),
        size: file.data_block.size,
        bucket_name: file.bucket.name,
        created_at: DateTime.now
      )
    rescue => e
      Metis.instance.logger.error("Failed to log unlink event: #{e.message}")
    end

    def self.log_vacuum(datablock:, project_name:, user:)
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

      # Orphaned = (linked by this project) ∩ (unlinked) - (currently in use) - (already vacuumed)
      orphaned_ids = (linked_datablock_ids & unlinked_datablock_ids) - used_datablock_ids - vacuumed_datablock_ids

      Metis::DataBlock
        .where(id: orphaned_ids)
        .exclude(removed: true)
        .all
    end
  end
end

