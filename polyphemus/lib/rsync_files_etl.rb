require "rsync"
require_relative "etl"
require_relative "filename_scan_based_etl_scanner"

class Polyphemus
  class RsyncFilesEtlCursor < EtlCursor
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

  # Abstract base class for an ETL that scans a remote
  #   server directory via rsync
  #   and synchronizes the list of available files with
  #   the Polyphemus database.
  class RsyncFilesEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_bucket_pairs:, host:, username:, password:, root:, scanner:, cursor_env:)
      @remote_path = "#{username}@#{host}:#{root}"

      cursors = cursors_from_pairs(
        pairs: project_bucket_pairs,
        pair_keys: %w[project_name bucket_name],
        cls: RsyncFilesEtlCursor,
        cursor_env: cursor_env
      )

      scanner = build_scanner if scanner.nil?

      super(
        cursors: cursors,
        scanner: scanner,
      )
    end

    def build_scanner
      FilenameScanBasedEtlScanner.new.result_updated_at do |change|
        change.timestamp
      end.result_id do |change|
        "#{host}://#{change.filename}"
      end.execute_batch_find do
        # RSync just lists everything at once, so
        #   we don't increment or batch any requests.
        results = Rsync.run(
          @remote_path,
          ".",
          ["-avzh",
            "--dry-run",
            "--itemize-changes",
            "--exclude=test",
            "--exclude=Reports",
            "--exclude=Stats",
            "--rsh=\"/usr/bin/sshpass -p #{password} ssh -l #{username}\""]
        )

        raise results.error unless results.success?

        results.changes.select do |change|
          :file == change.file_type
        end
      end
    end
  end
end
