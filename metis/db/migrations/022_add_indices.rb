Sequel.migration do
    change do
      alter_table(:files) do
        add_index :bucket_id
        add_index :folder_id
        add_index :data_block_id
      end

      alter_table(:folders) do
        add_index :bucket_id
        add_index :folder_id
      end
    end
  end
