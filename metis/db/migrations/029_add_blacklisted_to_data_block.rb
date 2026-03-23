Sequel.migration do
  up do
    alter_table(:data_blocks) do
      add_column :blacklisted, TrueClass, default: false
    end

    # Migrate removed: true -> blacklisted: true, removed: false
    from(:data_blocks).where(removed: true).update(blacklisted: true, removed: false)
  end

  down do
    # Reverse: blacklisted: true -> removed: true
    from(:data_blocks).where(blacklisted: true).update(removed: true)

    alter_table(:data_blocks) do
      drop_column :blacklisted
    end
  end
end
