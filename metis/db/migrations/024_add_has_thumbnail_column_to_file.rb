Sequel.migration do
  change do
    alter_table(:files) do
      add_column :has_thumbnail, TrueClass
    end
  end
end
