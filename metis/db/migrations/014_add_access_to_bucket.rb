Sequel.migration do
  change do
    alter_table(:buckets) do
      add_column :description, String
      add_column :owner, String
      add_column :access, String
    end

    self[:buckets].update(
      owner: 'metis',
      access: 'viewer'
    )

    alter_table(:buckets) do
      set_column_not_null :owner
      set_column_not_null :access
    end
  end
end
