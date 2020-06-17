Sequel.migration do 
  up do
    extension(:constraint_validations)
    create_constraint_validations_table

    create_table(:groups) do
      primary_key(:id)
      String(:group_name, {null: false})
    end

    create_table(:projects) do
      primary_key(:id)
      foreign_key(:group_id, :groups)
      String(:project_name, {null: false, unique: true})
      String(:project_name_full, {null: false, unique: true})
      String(:project_description)
    end

    create_table(:users) do
      primary_key(:id)
      String(:email, {null: false, unique: true})
      String(:first_name)
      String(:last_name)
      String(:pass_hash)
      DateTime(:user_create_stamp, {null: false, default: Time.now})
    end

    create_table(:apps) do
      primary_key(:id)
      String(:app_key, {null: false})
      String(:app_name, {null: false})
    end

    create_table(:permissions) do
      primary_key(:id)
      foreign_key(:user_id, :users)
      foreign_key(:project_id, :projects)

      String(:role, {null: false})
      validate do
        includes(['administrator', 'editor', 'viewer'], :role)
      end
    end

    create_table(:tokens) do
      primary_key(:id)
      String(:token, null: false)
      foreign_key(:user_id, :users)
      DateTime(:token_login_stamp, {null: false, default: Time.now})
      DateTime(:token_expire_stamp,{null: false})
      DateTime(:token_logout_stamp,{null: false})
    end
  end

  down do
    extension(:constraint_validations)
    drop_table(:tokens)
    drop_table(:permissions)
    drop_table(:apps)
    drop_table(:users)
    drop_table(:projects)
  end
end
