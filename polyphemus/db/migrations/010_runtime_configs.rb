Sequel.migration do
  up do
    create_table(:runtime_configs) do
      primary_key :id
      Integer :config_id, null: false, unique: true
      json :config, null: true
      Integer :run_interval, null: true
      Boolean :disabled, null: false, default: false
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    # Get all unique config_ids from configs table
    db = Polyphemus.instance.db
    config_ids = db[:configs].select(:config_id).distinct.all

    # Create runtime_config for each config_id
    config_ids.each do |config|
      db[:runtime_configs].insert(
        config_id: config[:config_id],
        created_at: Time.now,
        updated_at: Time.now
      )
    end
  end

  down do
    drop_table(:runtime_configs)
  end
end
