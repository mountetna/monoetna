Sequel.migration do
  change do
    alter_table(:data_blocks) do
      add_column :size, Integer
    end
  end
end
