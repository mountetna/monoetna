Sequel.migration do
  change do
    alter_table(:files) do
      drop_column :start_upload
      drop_column :finish_upload
      drop_column :hashing_algorithm
      rename_column :upload_by, :uploader
      rename_column :file_size, :size
    end
  end
end
