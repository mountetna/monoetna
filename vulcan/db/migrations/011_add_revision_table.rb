Sequel.migration do
  up do
    create_table(:revisions) do
      primary_key :id
      foreign_key :workspace_id, :workspace
      json :ui_params, null: false
      String :message, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    run('CREATE SEQUENCE revisions_ids START 1;')
  end

  down do
    drop_table(:revisions)

    run('DROP SEQUENCE revisions_ids')
  end
end