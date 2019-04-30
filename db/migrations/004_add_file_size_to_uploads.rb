Sequel.migration do
  change do
    alter_table(:uploads) do
      add_column :file_size, Integer, null: false
    end
  end
end
