Sequel.migration do
  up do
    alter_table(:users) do
      set_column_type :flags, :jsonb
    end
  end

  down do
    alter_table(:users) do
      set_column_type :flags, :json
    end
  end
end
