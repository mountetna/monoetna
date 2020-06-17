Sequel.migration do 
  change do
    alter_table(:users) do
      add_column :public_key, String
    end
  end
end
