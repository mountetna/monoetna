Sequel.migration do
  up do
    rename_table :etl_configs, :configs

    alter_table(:configs) do
      rename_column :etl, :workflow_name
      drop_column :params
      drop_column :ran_at
      drop_column :comment
      drop_column :status
      drop_column :output
      set_column_type :run_interval, String, size: 50
    end
  end

  down do
    alter_table(:configs) do
      rename_column :workflow_name, :etl
      add_column :configs, :params, :json, null: false
      add_column :configs, :ran_at, DateTime
      add_column :configs, :comment, String
      add_column :configs, :status, String
      add_column :configs, :output, String
      set_column_type :run_interval, Integer
    end

    rename_table :configs, :etl_configs
  end
end
