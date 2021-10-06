Sequel.migration do
  change do
    alter_table(:files) do
      add_index :updated_at
    end

    alter_table(:folders) do
      add_index :updated_at
    end
  end
end
