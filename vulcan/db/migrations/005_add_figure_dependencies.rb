Sequel.migration do
  change do
    alter_table(:figures) do
      add_column :dependencies, :json, null: false, default: "{}"
    end
  end
end
