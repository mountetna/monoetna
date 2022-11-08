Sequel.migration do
  change do
    alter_table(:identifiers) do
      drop_index [:identifier]
      add_index [:project_name, :identifier], unique: true
    end
  end
end
