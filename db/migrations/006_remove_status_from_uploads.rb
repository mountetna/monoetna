Sequel.migration do
  change do
    alter_table(:uploads) do
      drop_column :status
    end
  end
end
