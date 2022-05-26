Sequel.migration do
  change do
    create_table(:cc_agreements) do
      primary_key(:id)
      DateTime(:timestamp, { null: false, default: Time.now })
      String(:cc_text, { null: false })
      String(:project_name, { null: false })
      String(:user_email, { null: false })
      TrueClass(:agreed, { null: false })
    end
  end
end
