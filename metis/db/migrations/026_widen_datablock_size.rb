Sequel.migration do
  up do
    # for dev environments where the migration was not run by hand.
    unless Metis.instance.environment == :production
      from(:data_blocks).where(size: nil).each do |db|
        db.update(size: db.actual_size || 0)
      end
    end
    alter_table(:data_blocks) do
      set_column_type :size, :bigint, null: false
    end
  end
  down do
    alter_table(:data_blocks) do
      set_column_type :size, Integer, null: true
    end
  end
end
