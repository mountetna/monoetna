Sequel.migration do
  change do
    create_table(:workflow_errors) do
      String :hash, nil: false
      String :message, nil: false
      String :uuid, nil: false

      index [:hash], unique: true
    end
  end
end
