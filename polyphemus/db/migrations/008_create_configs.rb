Sequel.migration do
  up do
    create_table(:configs) do
      primary_key :id
      String :project_name, null: false
      String :workflow_name, null: false
      String :user, null: false
      String :name, null: false
      Integer :version_number, null: false
      column :config, :json, null: false
      String :secrets, null: false
      String :run_interval, size: 50
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:configs)
  end
end
