Sequel.migration do
  up do
    alter_table(:data_blocks) do
      add_column :removed, TrueClass, default: false
    end
  end

  down do
    alter_table(:data_blocks) do
      drop_column :removed
    end
  end
end
