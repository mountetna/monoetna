Sequel.migration do 
  up do
    create_table(:etls) do
      primary_key(:id)
      String :name, unique: true, null: false
      column :config, :json, null: false
      String :project_name, null: false
      DateTime :updated_at, null: false
      DateTime :created_at, null: false
      Integer :run_interval
    end
  end

  down do
    drop_table(:etls)
  end
end
