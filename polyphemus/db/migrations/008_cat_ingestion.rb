Sequel.migration do
  up do
    create_table(:cat_ingestion) do
      primary_key :id
      foreign_key :etl_configs_id, :etl_configs, null: false, on_delete: :cascade
      String :argo_id, null: false
      DateTime :last_scan, null: false
      Integer :num_files_to_update, null: false, default: 0
      Integer :num_c4_files_updated, null: false, default: 0
      Integer :num_metis_files_updates, null: false, default: 0
      DateTime :updated_at, null: false
      DateTime :modified_at, null: false
    end
  end

  down do
    drop_table(:cat_ingestion)
  end
end
