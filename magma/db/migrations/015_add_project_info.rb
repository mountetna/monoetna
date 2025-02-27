Sequel.migration do
  change do
    create_table(:enums) do
      primary_key :id
      String :type, null: false
      String :value, null: false
    end

    create_table(:project_infos) do
      primary_key :id
      String :project_name, null: false
      String :project_name_full, null: false
      String :principal_investigator, null: false
      String :email, null: false
      String :description, null: false
      String :funding_source, null: false
      TrueClass :data_collection_complete, null: false
      foreign_key :theme_id, :enums, null: false
      foreign_key :status_id, :enums, null: false
      DateTime :created_at
      DateTime :updated_at
    end

    create_table(:project_info_enums) do
      foreign_key :project_info_id
      foreign_key :enum_id
    end
  end
end
