Sequel.migration do
  up do
    alter_table(:uploads) do
      set_column_type :file_size, :bigint
      set_column_type :next_blob_size, :bigint
      set_column_type :current_byte_position, :bigint
    end
  end
  down do
    alter_table(:uploads) do
      set_column_type :file_size, Integer
      set_column_type :next_blob_size, Integer
      set_column_type :current_byte_position, Integer
    end
  end
end
