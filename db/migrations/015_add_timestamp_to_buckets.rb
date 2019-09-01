Sequel.migration do
  up do
    alter_table(:buckets) do
      add_column :created_at, DateTime
      add_column :updated_at, DateTime
    end

    from(:buckets).update(
      created_at: DateTime.now,
      updated_at: DateTime.now
    )

    alter_table(:buckets) do
      set_column_not_null :created_at
      set_column_not_null :updated_at
    end
  end

  down do
    alter_table(:buckets) do
      drop_column :created_at
      drop_column :updated_at
    end
  end
end
