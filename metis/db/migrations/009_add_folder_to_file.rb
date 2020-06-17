Sequel.migration do
  change do
    alter_table(:files) do
      add_foreign_key :folder_id, :files
      add_column :is_folder, TrueClass, default: false
    end
  end
end
