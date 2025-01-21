Sequel.migration do
  up do
    rename_table :etl_configs, :configs

    alter_table(:configs) do
      rename_column :etl, :workflow_type
      rename_column :name, :workflow_name
      drop_column :params 
      drop_column :ran_at
      drop_column :status
      drop_column :output
      drop_column :run_interval
      add_column :updated_at, DateTime
    end
  end

  down do
    alter_table(:configs) do
      rename_column :workflow_name, :etl
      add_column :params, :json, null: false
      add_column :ran_at, DateTime
      add_column :status, String
      add_column :output, String
      add_column :run_interval, Integer
      drop_column :updated_at
    end

    rename_table :configs, :etl_configs
  end
end
