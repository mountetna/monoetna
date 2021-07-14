require_relative "etl"
require_relative "time_scan_based_etl_scanner"

class Polyphemus
  class DbTriageFileEtlCursor < EtlCursor
    def initialize(job_name:, project_name:, bucket_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      super("#{job_name}_metis_triage_#{project_name}_#{bucket_name}")
      self[:project_name] = project_name
      self[:bucket_name] = bucket_name
    end

    def reset!
      super { self[:seen_ids] = [] }
    end
  end

  # Abstract base class for an ETL that scans the Polyphemus
  #   database for newly updated files that have been
  #   flagged for triage ingestion.
  class DbTriageFileEtl < Etl
    # Subclasses should provide default values here, since commands are constructed
    def initialize(project_bucket_pairs:, table_name: nil, limit: 20, timeout: nil)
      file_cursors = project_bucket_pairs.map do |project_name, bucket_name|
        DbTriageFileEtlCursor.new(job_name: self.class.name, project_name: project_name, bucket_name: bucket_name).load_from_db
      end

      @table_name = table_name
      @limit = limit
      @timeout = timeout

      super(
        cursor_group: EtlCursorGroup.new(file_cursors),
        scanner: TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
          Polyphemus.instance.db[@table_name.to_sym].where(
            should_ingest: true,
          )
        end.result_updated_at do |file|
          file[:updated_at]
        end.result_id do |file|
          "#{file[:host]}://#{file[:name]}"
        end.execute_batch_find do |query, i|
          query.limit(@limit * i)

          Polyphemus.instance.db.transaction do
            Polyphemus.instance.db.run("SET LOCAL statement_timeout = #{@timeout}") if @timeout
            Polyphemus.instance.db[query.sql].all
          end
        end,
      )
    end
  end
end
