Sequel.migration do
  change do
    alter_table(:users) do
      add_column :created_at, DateTime, null: false, default: Time.now
      add_column :updated_at, DateTime, null: false, default: Time.now
    end

    alter_table(:projects) do
      add_column :created_at, DateTime, null: false, default: Time.now
      add_column :updated_at, DateTime, null: false, default: Time.now
    end
  end
end
