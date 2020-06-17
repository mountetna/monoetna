Sequel.migration do
  change do
    alter_table(:files) do
      add_column :archive_id, String
    end
  end
end
