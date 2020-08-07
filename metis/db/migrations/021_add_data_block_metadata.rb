Sequel.migration do
  up do
    create_table(:data_block_metadatas) do
      primary_key :id

      foreign_key :data_block_id, :data_blocks, null: false
      String :key, null: false
      String :value, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    alter_table(:data_block_metadatas) do
      add_index [:key, :data_block_id], name: "data_block_metadata_key_data_block_id_idx", unique: true
    end
  end

  down do
    drop_table(:data_block_metadatas)
    drop_index "data_block_metadata_key_data_block_id_idx"
  end
end
