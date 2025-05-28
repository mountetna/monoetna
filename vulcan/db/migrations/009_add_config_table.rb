Sequel.migration do
  up do
    create_table(:configs) do
      primary_key :id
      foreign_key :workspace_id, :workspaces, null: false
      String :path, null: false
      String :hash, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:configs)
  end
end
