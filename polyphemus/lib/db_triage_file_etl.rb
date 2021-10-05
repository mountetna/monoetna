require_relative "db_triage_file_base_etl"

class Polyphemus
  # Base class for an ETL that scans the Polyphemus
  #   database for newly updated files that have been
  #   flagged for triage ingestion.
  class DbTriageFileEtl < DbTriageFileBaseEtl
    def initialize_query(cursor)
      Polyphemus::IngestFile.where(
        should_ingest: true,
        @column_name => nil,
        removed_from_source: false,
      )
    end
  end
end
