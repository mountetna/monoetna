Sequel.migration do
  change do
    alter_table(:uploads) do
      add_column :metis_uid, String, null: false
    end
  end
end
