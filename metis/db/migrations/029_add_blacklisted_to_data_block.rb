Sequel.migration do
  up do
    alter_table(:data_blocks) do
      rename_column :removed, :blacklisted
    end
  end

  down do
    alter_table(:data_blocks) do
      rename_column :blacklisted, :removed
    end
  end
end
