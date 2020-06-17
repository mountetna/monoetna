Sequel.migration do
  change do
    create_table(:folders) do
      primary_key :id
      String :project_name, text: true, null: false
      String :folder_name, text: true, null: false

      foreign_key :bucket_id, :buckets, key: [:id]
      foreign_key :folder_id, :folders, key: [:id]

      TrueClass :read_only, default: false, null: false
      String :author, text: true, null: false

      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
    alter_table(:files) do
      drop_column :folder_id
      drop_column :is_folder
      add_foreign_key :folder_id, :folders
    end
  end
end
