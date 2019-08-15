Sequel.migration do
  change do
    alter_table(:buckets) do
      add_column :description, String
      add_column :owner, String, null: false
      add_column :access, String, null: false
    end
  end
end
