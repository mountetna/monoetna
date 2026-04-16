Sequel.migration do
  up do
    alter_table(:models) do
      add_column :template_project_name, String, null: true
      add_column :template_model_name, String, null: true
      add_index [:template_project_name, :template_model_name]
      add_foreign_key [:template_project_name, :template_model_name], :models,
        key: [:project_name, :model_name],
        name: :models_template_reference_fkey
    end
  end

  down do
    alter_table(:models) do
      drop_constraint :models_template_reference_fkey, type: :foreign_key
      drop_index [:template_project_name, :template_model_name]
      drop_column :template_model_name
      drop_column :template_project_name
    end
  end
end
