Sequel.migration do
  up do
    create_table(:watch_folders) do
      primary_key(:id)
      String :project_name, null: false
      String :bucket_name, null: false
      String :folder_path, null: false
      Integer :metis_id, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
      String :watch_type, null: false
    end

    alter_table(:watch_folders) do
      add_index [:project_name, :bucket_name, :metis_id, :watch_type], :unique => true
    end
  end

  down do
    drop_table(:watch_folders)
  end
end
