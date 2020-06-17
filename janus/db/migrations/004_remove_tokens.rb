Sequel.migration do 
  change do
    drop_table(:tokens)
  end
end
