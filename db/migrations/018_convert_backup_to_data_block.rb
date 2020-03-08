Sequel.migration do
  up do
    rename_table(:backups, :data_blocks)
    alter_table(:data_blocks) do
      set_column_allow_null :archive_id
      add_column :created_at, DateTime
      add_column :updated_at, DateTime
      drop_constraint(:backups_description_key, type: :unique)
    end

    from(:data_blocks).update(
      created_at: DateTime.now,
      updated_at: DateTime.now
    )

    alter_table(:data_blocks) do
      set_column_not_null :created_at
      set_column_not_null :updated_at
    end

    alter_table(:files) do
      rename_column :backup_id, :data_block_id
      set_column_not_null :data_block_id
      drop_column :file_hash
    end
  end

  down do
    rename_table(:data_blocks, :backups)
    alter_table(:backups) do
      set_column_not_null :archive_id
      drop_column :created_at
      drop_column :updated_at
      add_unique_constraint :description
    end
    alter_table(:files) do
      rename_column :data_block_id, :backup_id
      add_column :file_hash, String
      set_column_allow_null :backup_id
    end
  end
end
