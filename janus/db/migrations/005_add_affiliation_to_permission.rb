Sequel.migration do 
  change do
    alter_table(:permissions) do
      add_column :affiliation, String
    end
  end
end
