Sequel.migration do
  change do
    create_table(:flags) do
      primary_key :id
      String :project_name, null: false
      String :flag_name, null: false
      column :possible_values, "text[]", null: false
      String :value, null: true
    end
  end
end