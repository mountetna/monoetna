Sequel.migration do 
  up do
    create_table(:etl_configs) do
      primary_key(:id)
      String :project_name, null: false
      String :name, null: false
      String :etl, null: false
      column :config, :json, null: false
      String :secrets, null: false
      DateTime :updated_at, null: false
      DateTime :created_at, null: false
      DateTime :ran_at
      String :comment
      String :status
      String :output
      Integer :run_interval, null: false
      TrueClass :archived, null: false, default: false
    end
  end

  down do
    drop_table(:etl_configs)
  end
end
