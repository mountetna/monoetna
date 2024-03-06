Sequel.migration do
  up do
    create_table(:runs) do
      primary_key :id
      Integer :workspace_id, null: false
      foreign_key [:workspace_id], :workspaces
      column :outputs_generated, 'text[]', null: false
      column :jobs_invoked, 'text[]', null: false
    end

    run('CREATE SEQUENCE runs_ids START 1;')
  end

  down do
    drop_table(:runs)

    run('DROP SEQUENCE runs_ids')
  end
end
