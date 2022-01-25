require_relative "db_triage_file_base_etl"

class Polyphemus
  # Base class for an ETL that scans the Polyphemus
  #   database for newly updated files that have been
  #   ingested already.
  class DbTriageFileNotificationEtl < DbTriageFileBaseEtl
    # Batch notifications in groups of 100 files
    def initialize(
      project_bucket_pairs:,
      column_name:,
      cursor_env: {},
      scanner: build_scanner,
      limit: 100,
      timeout: nil
    )
      super
    end

    def initialize_query(cursor)
      # Careful, apparently exclude(condA, condB) => exclude where BOTH condA and condB,
      # but we want exclude where EITHER condA and condB, so we're forced to make two separate calls.
      Polyphemus::IngestFile.exclude(
        @column_name => nil,
      ).exclude(
        :removed_from_source => true,
      ).where { updated_at > cursor.updated_at }
    end
  end
end
