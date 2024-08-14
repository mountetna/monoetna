Sequel.migration do
  up do
    create_table(:runs) do
      primary_key :id
      Integer :workspace_id, null: false
      foreign_key [:workspace_id], :workspaces
      String :slurm_run_uuid, null: false
      String :log_path, null: false
      String :config_path, null: false
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
