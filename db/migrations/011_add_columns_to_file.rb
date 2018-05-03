Sequel.migration do
  up do
    alter_table(:files) do
      add_column :author, String
      add_column :created_at, DateTime
      add_column :updated_at, DateTime
    end

    from(:files).update(
      author: '',
      created_at: DateTime.now,
      updated_at: DateTime.now
    )

    alter_table(:files) do
      set_column_not_null :author
      set_column_not_null :created_at
      set_column_not_null :updated_at
    end

    alter_table(:uploads) do
      # we don't explicitly link to the file any more
      drop_column :file_id

      # add the file attributes here
      add_column :file_name, String, null: false
      add_column :project_name, String, null: false
      add_column :author, String, null: false
      add_foreign_key :bucket_id, :buckets
    end
  end

  down do
    alter_table(:files) do
      drop_column :author
      drop_column :created_at
      drop_column :updated_at
    end
    alter_table(:uploads) do
      drop_column :author
      drop_column :file_name
      drop_column :project_name
      drop_column :bucket_id
    end
  end
end
