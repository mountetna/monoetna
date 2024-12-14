Sequel.migration do
  up do
    create_table(:runs) do
      primary_key :id
      String :run_id, null: false
      json :state, null: false
      Integer :config_id, null: false
      Integer :version_number, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
    # our version of Sequel does not support foreign keys on composite keys
    run(
      "ALTER TABLE runs " \
      "ADD CONSTRAINT runs_config_version_number_fkey " \
      "FOREIGN KEY (config_id, version_number) " \
      "REFERENCES configs (config_id, version_number);"
    )
  end

  down do
    run("ALTER TABLE runs DROP CONSTRAINT runs_config_version_number_fkey;")
    drop_table(:runs)
  end
end
