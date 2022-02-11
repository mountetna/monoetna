require_relative "etl"
require_relative "time_scan_based_etl_scanner"

class Polyphemus
  class DbTriageFileBaseEtlCursor < EtlCursor
    def initialize(job_name:, project_name:, bucket_name:)
      raise "project_name cannot be nil" if project_name.nil?
      raise "bucket_name cannot be nil" if bucket_name.nil?
      raise "batch_end_at must be set if updated_at is set." if updated_at && !batch_end_at
      raise "updated_at must be set if batch_end_at is set." if batch_end_at && !updated_at
      super("#{job_name}_triage_ingest_#{project_name}_#{bucket_name}")
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
  class DbTriageFileBaseEtl < Etl
    def initialize(project_bucket_pairs:, column_name:, cursor_env: {}, scanner: build_scanner, limit: 20, timeout: nil)
      @limit = limit
      @timeout = timeout
      @column_name = column_name
      @project_bucket_pairs = project_bucket_pairs

      cursors = cursors_from_pairs(
        pairs: project_bucket_pairs,
        pair_keys: %w[project_name bucket_name],
        cls: DbTriageFileBaseEtlCursor,
        cursor_env: cursor_env
      )

      super(cursors: cursors, scanner: scanner)
    end

    def initialize_query(cursor)
      raise "subclasses must implement initialize_query"
    end

    def build_scanner
      TimeScanBasedEtlScanner.new.start_batch_state do |cursor|
        query = initialize_query(cursor)
        if (end_at = cursor[:batch_end_at])
          query = query.where {
            updated_at <= end_at + 1
          }
        end

        query
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
      end
    end
  end
end
