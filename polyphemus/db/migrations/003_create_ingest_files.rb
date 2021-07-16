Sequel.migration do
  up do
    create_table(:ingest_files) do
      primary_key(:id)
      String :name, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
      DateTime :ingested_at
      Boolean :should_ingest, null: false
      String :host, null: false
    end

    alter_table(:ingest_files) do
      add_index [:name, :host], :unique => true
    end
  end

  down do
    drop_table(:ingest_files)
  end
end
