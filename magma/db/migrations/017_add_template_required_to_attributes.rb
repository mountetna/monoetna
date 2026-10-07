Sequel.migration do
  up do
    alter_table(:attributes) do
      add_column :template_required, TrueClass, null: false, default: false
    end
  end

  down do
    alter_table(:attributes) do
      drop_column :template_required
    end
  end
end
