Sequel.migration do
  change do
    alter_table(:files) do
      rename_column :hash, :file_hash
    end
  end
end
