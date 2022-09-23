Sequel.migration do
  change do
    create_table(:rules) do
      primary_key :id

      String :project_name, null: false
      Integer :version_number, null: false
      column :config, :json, null: false
      DateTime :created_at, null: false

      index [:project_name, :version_number], unique: true
    end

    create_table(:identifiers) do
      primary_key :id

      String :project_name, null: false
      String :identifier, null: false
      String :author, null: false
      foreign_key :rule_id, :rules, null: false
      foreign_key :renamed_to_id, :identifiers
      DateTime :created_at, null: false

      index [:identifier], unique: true
    end
  end
end
