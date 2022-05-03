Sequel.migration do
  change do
    alter_table(:figures) do
      add_column :archived, TrueClass, null: false, default: false
      add_column :comment, String, null: false, default: 'Initial commit'
    end
  end
end
