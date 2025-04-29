Sequel.migration do
  change do
    create_table(:query_histories) do
      primary_key :id
      String :project_name, null: false
      String :user, null: false
      String :comment, null: false
      String :query, null: false
      DateTime :created_at, null: false
    end
  end
end
