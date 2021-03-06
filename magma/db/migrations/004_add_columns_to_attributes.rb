Sequel.migration do
  up do
    alter_table(:attributes) do
      add_column :type, String, null: false
      add_column :hidden, "boolean", default: false
      add_column :read_only, "boolean", default: false
      add_column :unique, "boolean", default: false
      add_column :index, "boolean", default: false
      add_column :restricted, "boolean", default: false

      set_column_not_null :project_name
      set_column_not_null :model_name
      set_column_not_null :attribute_name
    end
  end

  down do
    alter_table(:attributes) do
      drop_column :type
      drop_column :hidden
      drop_column :read_only
      drop_column :unique
      drop_column :index
      drop_column :restricted

      set_column_allow_null :project_name
      set_column_allow_null :model_name
      set_column_allow_null :attribute_name
    end
  end
end
