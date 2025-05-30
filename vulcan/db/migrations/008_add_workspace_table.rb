Sequel.migration do
  up do
    create_table(:workspaces) do
      primary_key :id

      Integer :workflow_id, null: false
      foreign_key [:workflow_id], :workflows
      json :target_mapping, null: false
      String :name, null: false
      String :user_email, null: false
      String :path, null: false
      String :git_version, null: false
      column :tags, 'text[]'
      column :dag, 'text[]'

      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    run('CREATE SEQUENCE workspaces_ids START 1;')
  end

  down do
    drop_table(:workspaces)

    run('DROP SEQUENCE workspaces_ids')
  end
end
