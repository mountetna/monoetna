Sequel.migration do
  change do
    alter_table(:projects) do
      add_column :requires_agreement, TrueClass, default: false
      add_column :cc_text, String, null: false, default: ""
      add_column :contact_email, String, null: false, default: ""
    end
  end
end
