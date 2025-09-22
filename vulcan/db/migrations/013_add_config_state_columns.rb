Sequel.migration do
  up do
    alter_table(:configs) do
      add_column :input_files, 'text[]'
      add_column :input_params, :jsonb
      add_column :future_state, :jsonb
    end
  end

  down do
    alter_table(:configs) do
      drop_column :input_files
      drop_column :input_params
      drop_column :future_state
    end
  end
end
