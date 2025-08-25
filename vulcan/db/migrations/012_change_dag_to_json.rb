Sequel.migration do
  up do
    alter_table(:workspaces) do
      # Change dag column from text[] to json
      drop_column :dag
      add_column :dag, :json
    end
  end

  down do
    alter_table(:workspaces) do
      # Revert back to text[] if needed
      drop_column :dag
      add_column :dag, 'text[]'
    end
  end
end 