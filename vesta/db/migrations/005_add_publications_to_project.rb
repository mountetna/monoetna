Sequel.migration do
  change do
    alter_table(:projects) do
      add_column :publications, :json, null: false, default: Sequel.pg_json_wrap([])
    end
  end
end
