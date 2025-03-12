Sequel.migration do
  up do
    create_table(:logs) do
      primary_key :id
      String :application, null: false
      String :project_name, null: false
      String :user, null: true
      String :event, null: false
      String :message, null: false
      json :payload, null: true
      DateTime :created_at, null: false
      Boolean :hidden, null: false, default: false
    end
  end

  down do
    drop_table(:logs)
  end
end
