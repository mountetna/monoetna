Sequel.migration do
  change do
    alter_table(:permissions) do
      add_index [:user_id, :project_id], unique: true
    end
  end
end
