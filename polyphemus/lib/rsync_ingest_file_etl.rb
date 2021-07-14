require "rsync"
require_relative "etl"
require_relative "memoryless_etl_scanner"

class Polyphemus
  class RsyncIngestFileEtlCursor < EtlCursor
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
  class RsyncIngestFileEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_bucket_pairs:, host:, username:, password:, root:)
      file_cursors = project_bucket_pairs.map do |project_name, bucket_name|
        RsyncIngestFileEtlCursor.new(job_name: self.class.name, project_name: project_name, bucket_name: bucket_name).load_from_db
      end

      remote_path = "#{username}@#{host}:#{root}"

      super(
        cursor_group: EtlCursorGroup.new(file_cursors),
        scanner: MemorylessEtlScanner.new.execute_batch_find do
          # RSync just lists everything at once, so
          #   we don't increment or batch any requests.
          binding.pry
          results = Rsync.run(
            remote_path,
            ".",
            ["-avzh",
             "--dry-run",
             "--itemize-changes",
             "--exclude=test",
             "--exclude=Reports",
             "--exclude=Stats",
             "--rsh=\"/usr/bin/sshpass -p #{password} ssh -o StrictHostKeyChecking=no -l #{username}\""]
          )

          raise results.error unless results.success?

          results.changes
        end,
      )
    end
  end
end
