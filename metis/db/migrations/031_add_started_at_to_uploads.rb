Sequel.migration do
  change do
    alter_table(:uploads) do
      add_column :started_at, DateTime
    end
  end
end
