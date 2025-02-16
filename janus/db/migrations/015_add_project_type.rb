Sequel.migration do
  up do
    alter_table(:projects) do
      add_column :project_type, String
    end

    db = Janus.instance.db

    db[:projects].where(requires_agreement: true, resource: true).update(
      project_type: "community"
    )
    db[:projects].where(requires_agreement: false, resource: true).update(
      project_type: "resource"
    )
    db[:projects].where(resource: false).update(
      project_type: "team"
    )

    alter_table(:projects) do
      drop_column :requires_agreement
      drop_column :resource
      set_column_not_null :project_type
      set_column_default :project_type, 'team'
    end
  end

  down do
    alter_table(:projects) do
      add_column :resource, TrueClass, default: false
      add_column :requires_agreement, TrueClass, default: false
    end

    db = Janus.instance.db

    db[:projects].where(project_type: 'community').update(
      resource: true,
      requires_agreement: true
    )
    db[:projects].where(project_type: 'resource').update(
      resource: true,
      requires_agreement: false
    )
    db[:projects].where(project_type: 'team').update(
      resource: false,
      requires_agreement: false
    )

    alter_table(:projects) do
      drop_column :project_type
    end
  end
end
