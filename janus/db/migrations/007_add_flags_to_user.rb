Sequel.migration do 
  change do
    alter_table(:users) do
      add_column :flags, :json
    end
  end
end
