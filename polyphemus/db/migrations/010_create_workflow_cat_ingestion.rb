Sequel.migration do
  up do
    create_table(:workflow_cat_ingestion) do
      primary_key :id
      String :run_id, null: false
      foreign_key :configs_id, :configs, null: false, on_delete: :cascade
      column :state, :json, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:cat_ingestion)
  end
end
