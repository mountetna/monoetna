Sequel.migration do
  change do
    create_table(:data_block_ledger) do
      primary_key :id
      String :project_name, null: true
      String :md5_hash, null: false, size: 40 
      String :file_path, text: true
      Integer :file_id
      Integer :data_block_id, null: false
      String :event_type, null: false, size: 40
      DateTime :created_at, null: false
      String :triggered_by, text: true
      column :size, 'bigint'
      String :bucket_name
      
      index [:data_block_id], name: :data_block_ledger_block_id
      index [:project_name, :event_type], name: :data_block_ledger_project_event
      index [:triggered_by, :event_type], name: :data_block_ledger_trigger_event
    end
  end
end

