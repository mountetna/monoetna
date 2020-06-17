Sequel.migration do 
  change do
    alter_table(:permissions) do
      add_column :privileged, TrueClass
    end
  end
end
