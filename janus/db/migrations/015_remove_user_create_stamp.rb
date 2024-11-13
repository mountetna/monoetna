Sequel.migration do
  up do
    alter_table(:users) do
      drop_column :user_create_stamp
    end
  end

  down do
    alter_table(:users) do
      add_column :user_create_stamp, DateTime, null: false, default: Sequel::CURRENT_TIMESTAMP
    end
  end
end
