Sequel.migration do
  up do
    create_table(:workflow_run_metadata) do
      primary_key :id
      String :run_id, null: false
      foreign_key :configs_id, :configs, null: false, on_delete: :cascade
      String :config_id, null: false
      column :meta_data, :json, null: false
      Text :comment, text: true
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:workflow_run_metadata)
  end
end
