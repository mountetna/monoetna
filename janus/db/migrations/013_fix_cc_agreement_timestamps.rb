Sequel.migration do
  change do
    alter_table(:cc_agreements) do
      add_column :created_at, DateTime
      add_column :updated_at, DateTime
      drop_column :timestamp
    end
  end
end
