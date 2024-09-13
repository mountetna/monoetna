Sequel.migration do
  up do
    create_table(:workflows) do
      primary_key :id
      column :projects, 'text[]', null: false
      String :name, null: false
      String :branch, null: false
      String :repo_remote_url, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    run('CREATE SEQUENCE workflows_ids START 1;')
  end

  down do
    drop_table(:workflows)

    run('DROP SEQUENCE workflows_ids')
  end
end
