Sequel.migration do
  up do
    alter_table(:runs) do
      add_column :input_files, 'text[]'
    end
  end

  down do
    alter_table(:runs) do
      drop_column :input_files
    end
  end
end
