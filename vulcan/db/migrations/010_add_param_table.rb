Sequel.migration do
  up do
    create_table(:params) do
      primary_key :id
      foreign_key :run_id, :runs
      json :data, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    run('CREATE SEQUENCE params_ids START 1;')
  end

  down do
    drop_table(:params)

    run('DROP SEQUENCE params_ids')
  end
end
