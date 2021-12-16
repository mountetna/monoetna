Sequel.migration do
  change do
    alter_table(:data_blocks) do
      add_column :has_thumbnail, TrueClass
    end
  end
end
