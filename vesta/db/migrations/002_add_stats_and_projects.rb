Sequel.migration do
  up do
    create_table(:projects) do
      primary_key :id

      String :name, null: false, index: true
      String :full_name, null: false
      # not adding indexes here since
      # this table is temporary and will be
      # replaced by direct magma queries
      String :description, null: true
      String :funding_source, null: true
      String :status, null: false
      String :type, null: false
      String :species, null: false
      DateTime :start_date, null: false
      String :theme, null: false
      Boolean :data_collection_complete, null: false
      DateTime :created_at, null: false, default: Sequel::CURRENT_TIMESTAMP
      DateTime :updated_at, null: false, default: Sequel::CURRENT_TIMESTAMP
    end

    # TODO: figure out why these columns can't be added as above
    alter_table(:projects) do
      add_column(:principal_investigators, JSON, null: false)
      add_column(:data_types, JSON, null: false)
    end

    create_table(:project_stats) do
      primary_key :id

      String :name, null: false, index: true
      Integer :user_count, null: false
      Integer :file_count, null: false
      Integer :byte_count, null: false
      Integer :subject_count, null: false
      Integer :sample_count, null: false
      Integer :assay_count, null: false
      Integer :clinical_data_count, null: false
      DateTime :recorded_at, null: false, index: true, default: Sequel::CURRENT_TIMESTAMP
    end

    create_table(:global_stats) do
      primary_key :id

      Integer :user_count, null: false
      Integer :file_count, null: false
      Integer :byte_count, null: false
      Integer :subject_count, null: false
      Integer :sample_count, null: false
      Integer :assay_count, null: false
      DateTime :recorded_at, null: false, index: true, default: Sequel::CURRENT_TIMESTAMP
    end
  end

  down do
    drop_table(:global_stats)
    drop_table(:project_stats)
    drop_table(:projects)
  end
end
