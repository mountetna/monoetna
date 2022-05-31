Sequel.migration do
  up do
    extension(:constraint_validations)
    alter_table(:permissions) do
      drop_constraint("permissions_role_check")
      validate do
        includes(["administrator", "editor", "viewer", "guest"], :role, name: "permission_allowable_roles")
      end
    end
  end

  down do
    extension(:constraint_validations)
    alter_table(:permissions) do
      validate do
        drop_constraint("permission_allowable_roles")
        includes(["administrator", "editor", "viewer"], :role, name: "permissions_role_check")
      end
    end
  end
end
