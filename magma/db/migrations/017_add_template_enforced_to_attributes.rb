Sequel.migration do
  up do
    alter_table(:attributes) do
      add_column :template_enforced, TrueClass, null: false, default: false
    end
  end

  down do
    alter_table(:attributes) do
      drop_column :template_enforced
    end
  end
end
