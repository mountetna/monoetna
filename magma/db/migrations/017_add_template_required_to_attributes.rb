Sequel.migration do
  up do
    alter_table(:attributes) do
      add_column :required_on_child_models, TrueClass, null: false, default: false
    end
  end

  down do
    alter_table(:attributes) do
      drop_column :required_on_child_models
    end
  end
end
