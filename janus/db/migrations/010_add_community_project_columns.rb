Sequel.migration do
  change do
    alter_table(:projects) do
      add_column :requires_agreement, TrueClass, default: false
      add_column :cc_text, String, null: false, default: ''
      add_column :contact_email, String, null: false, default: ''
    end

    create_table(:cc_agreement) do
      primary_key(:id)
      DateTime(:timestamp, {null: false, default: Time.now})
      String(:cc_text, {null: false})
      String(:project_name, {null: false})
      String(:user_email, {null: false})
      TrueClass(:agreed, {null: false})
    end
  end
end
