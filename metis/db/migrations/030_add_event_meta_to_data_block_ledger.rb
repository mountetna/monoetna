Sequel.migration do
  change do
    alter_table(:data_block_ledger) do
      add_column :event_meta, :jsonb
    end
  end
end
