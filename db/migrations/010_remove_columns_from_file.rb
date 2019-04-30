Sequel.migration do
  change do
    alter_table(:files) do
      drop_column :size
      drop_column :uploader
      drop_column :original_name
      add_column :read_only, TrueClass, null: false, default: false
    end
  end
end
