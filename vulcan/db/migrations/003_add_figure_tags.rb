Sequel.migration do
  change do
    alter_table(:figures) do
      add_column :tags, :json
    end
  end
end
