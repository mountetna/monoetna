Sequel.migration do
  up do
    rename_table :etl_configs, :configs

    alter_table(:configs) do
      rename_column :etl, :workflow_name
      drop_columns :params, :ran_at, :comment, :status, :output
      set_column_type :run_interval, String, size: 50
    end
  end

  down do
    alter_table(:configs) do
      rename_column :workflow_name, :etl
      add_column :params, :json, null: false
      add_column :ran_at, DateTime
      add_column :comment, String
      add_column :status, String
      add_column :output, String
      set_column_type :run_interval, Integer
    end

    rename_table :configs, :etl_configs
  end
end
