require "rsync"
require_relative "etl"
require_relative "time_scan_based_etl_scanner"

class Polyphemus
  class SftpIngestFileEtlCursor < EtlCursor
    def initialize(job_name:, project_name:, bucket_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      super("#{job_name}_sftp_sync_#{project_name}_#{bucket_name}")
      self[:project_name] = project_name
      self[:bucket_name] = bucket_name
    end

    def reset!
      super { self[:seen_ids] = [] }
    end
  end

  # Abstract base class for an ETL that scans an SFTP server
  #   and synchronizes the list of available files with
  #   the Polyphemus database.
  class SftpIngestFileEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_bucket_pairs:, host:, username:, root:, limit: 20, timeout: nil)
      file_cursors = project_bucket_pairs.map do |project_name, bucket_name|
        SftpIngestFileEtlCursor.new(job_name: self.class.name, project_name: project_name, bucket_name: bucket_name).load_from_db
      end

      @limit = limit
      @timeout = timeout
      remote_path = "#{username}@#{host}:#{root}"

      super(
        cursor_group: EtlCursorGroup.new(file_cursors),
        scanner: TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
        end.result_updated_at do |change|
          change.timestamp
        end.result_id do |change|
          change.filename
        end.execute_batch_find do |query, i|
          # RSync just lists everything at once, so
          #   we don't increment or batch any requests.
          Rsync.run(
            remote_path,
            ["-avzhi",
             "--list-only",
             "--exclude=test",
             "--exclude=Reports",
             "--exclude=Stats"]
          )
        end,
      )
    end
  end
end
