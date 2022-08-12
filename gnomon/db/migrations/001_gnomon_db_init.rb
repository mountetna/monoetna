Sequel.migration do 
  up do
    create_table(:configs) do
      primary_key(:id)
      String :project_name, unique: true, null: false
      column :config, :json, null: false
      Integer :version, null: false, default: 0
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
  end
end
