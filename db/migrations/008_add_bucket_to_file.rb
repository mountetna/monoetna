Sequel.migration do
  change do
    create_table(:buckets) do
      primary_key :id
      String :name, null: false
      String :project_name, null: false
    end
    alter_table(:files) do
      add_foreign_key :bucket_id, :buckets
    end
  end
end
