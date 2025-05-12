Sequel.migration do
  change do
    alter_table(:files) do
      add_column :restricted, TrueClass
    end
  end
end
