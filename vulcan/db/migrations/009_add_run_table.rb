Sequel.migration do
  up do
    create_table(:runs) do
      primary_key :id
      Integer :workspace_id, null: false
      foreign_key [:workspace_id], :workspaces
      Json :run_config, null: false
      Json :job_status, null: true
      String :executor, null: false
      Integer :executor_job_id, null: false
      column :outputs_generated, 'text[]', null: true
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    run('CREATE SEQUENCE runs_ids START 1;')
  end

  down do
    drop_table(:runs)

    run('DROP SEQUENCE runs_ids')
  end
end
