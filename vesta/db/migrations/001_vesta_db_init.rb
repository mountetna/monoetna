Sequel.migration do
  up do
    extension(:constraint_validations)
    create_constraint_validations_table
    extension(:pg_array)
    extension(:pg_json)
  end

  down do
    extension(:pg_json)
    extension(:pg_array)
    extension(:constraint_validations)
  end
end
