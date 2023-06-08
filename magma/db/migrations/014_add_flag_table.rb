Sequel.migration do
  change do
    create_table(:flags) do
      primary_key :id
      String :project_name, null: false
      String :flag_name, null: false
      String :value, null: false
      DateTime :created_at
      DateTime :updated_at
    end
  end
end