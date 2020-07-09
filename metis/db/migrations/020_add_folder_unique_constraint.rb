Sequel.migration do
  up do
    alter_table(:folders) do
      add_index [:bucket_id, :folder_name], name: "unique_folder_name_per_parent_null", where: 'folder_id is NULL', unique: true
      add_unique_constraint [:bucket_id, :folder_id, :folder_name], unique_constraint_name: "unique_folder_name_per_parent"
    end

  end

  down do
    alter_table(:folders) do
      drop_constraint("unique_folder_name_per_parent_null", :type=>:unique)
    end
  end
end
