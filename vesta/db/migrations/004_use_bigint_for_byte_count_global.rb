Sequel.migration do
  up do
    alter_table(:global_stats) do
      set_column_type :byte_count, :bigint
    end
  end
  down do
    alter_table(:global_stats) do
      set_column_type :byte_count, Integer
    end
  end
end