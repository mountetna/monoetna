Sequel.migration do
    change do
      alter_table(:files) do
        add_index :bucket_id
        add_index :folder_id
        add_index :data_block_id
        add_index :project_name
      end

      alter_table(:folders) do
        add_index :bucket_id
        add_index :folder_id
        add_index :project_name
      end

      alter_table(:buckets) do
        add_index :project_name
      end
    end
  end
