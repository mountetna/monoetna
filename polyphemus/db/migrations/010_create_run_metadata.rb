Sequel.migration do
  up do
    create_table(:run_metadata) do
      primary_key :id
      String :run_id, null: false
      json :orchestrator_metadata, null: true
      json :runtime_config, null: true
      String :output, null: true
      Integer :run_interval, null: false
      Integer :config_id, null: false
      Integer :version_number, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
    run(
      "ALTER TABLE run_metadata " \
      "ADD CONSTRAINT run_metadata_config_version_number_fkey " \
      "FOREIGN KEY (config_id, version_number) " \
      "REFERENCES configs (config_id, version_number);"
    )
  end

  down do
    run("ALTER TABLE run_metadata DROP CONSTRAINT run_metadata_config_version_number_fkey;")
    drop_table(:run_metadata)
  end
end
