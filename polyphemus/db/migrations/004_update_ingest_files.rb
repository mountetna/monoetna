Sequel.migration do
  change do
    alter_table(:ingest_files) do
      rename_column :ingested_at, :archive_ingested_at
      add_column :triage_ingested_at, DateTime
    end
  end
end
