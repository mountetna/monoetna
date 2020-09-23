Sequel.migration do 
  up do
    create_table(:cursors) do
      primary_key(:id)
      String :name, unique: true, null: false
      column :value, :json, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:cursors)
  end
end
