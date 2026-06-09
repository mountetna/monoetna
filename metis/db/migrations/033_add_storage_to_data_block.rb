Sequel.migration do
  up do
    alter_table(:data_blocks) do
      add_column :storage, String, default: 'default'
    end
  end

  down do
    alter_table(:data_blocks) do
      drop_column :storage
    end
  end
end
