Sequel.migration do
  up do
    create_table(:runs) do
      primary_key :id
      String :run_id, null: false
      foreign_key :config_id, :configs, null: false
      json :state, null: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:runs)
  end
end
