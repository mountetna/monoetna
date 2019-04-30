Sequel.migration do
  change do
    alter_table(:uploads) do
      drop_column :current_blob_size
    end
  end
end
