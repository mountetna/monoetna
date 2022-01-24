Sequel.migration do
  change do
    create_table(:sessions) do
      String :owner, nil: false
      json :params
      String :key, nil: false

      index [:key], unique: true
    end
  end
end
