Sequel.migration do
  up do
    create_table(:runtime_configs) do
      primary_key :id
      Integer :config_id, null: false, unique: true
      json :config, null: true
      Integer :run_interval, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end
  # TO think about: maybe foreign key to configs.config_id

  down do
    drop_table(:runtime_configs)
  end
end
