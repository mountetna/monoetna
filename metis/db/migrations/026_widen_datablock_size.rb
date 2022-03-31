Sequel.migration do
  up do
    alter_table(:data_blocks) do
      set_column_type :size, :bigint
    end
  end
  down do
    alter_table(:data_blocks) do
      set_column_type :size, Integer
    end
  end
end
