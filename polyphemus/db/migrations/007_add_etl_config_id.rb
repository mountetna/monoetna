Sequel.migration do 
  up do
    alter_table(:etl_configs) do
      add_column :config_id, Integer
      add_column :version_number, Integer
      drop_column :archived

      add_constraint(nil) {version_number > 0}

      add_unique_constraint([:config_id, :version_number])
    end

    run('CREATE SEQUENCE etl_configs_ids START 1;')

    # backfill ids for existing names

    db = Polyphemus.instance.db

    db[:etl_configs].all.group_by { |e| [ e[:project_name], e[:name] ] }.each do |(project_name,name), configs|
      configs = configs.sort_by{ |e| e[:updated_at ] }

      config_id = db.get { nextval("etl_configs_ids") }

      configs.each.with_index do |config,i|
        db[:etl_configs].where(
          id: config[:id]
        ).update(
          config_id: config_id,
          version_number: i+1
        )
      end
    end

    alter_table(:etl_configs) do
      set_column_not_null :config_id
      set_column_not_null :version_number
      drop_column :updated_at
    end
  end

  down do
    alter_table(:etl_configs) do
      drop_column :config_id
      drop_column :version_number
      add_column :updated_at, DateTime, null: false
      add_column :archived, TrueClass, null: false, default: false
    end

    run('DROP SEQUENCE etl_configs_ids')
  end
end
