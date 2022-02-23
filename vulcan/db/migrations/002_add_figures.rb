Sequel.migration do
  up do
    create_table(:figures) do
      primary_key :id

      Integer :figure_id, null: false

      String :author, null: false
      String :title, null: false
      String :project_name, null: false
      String :workflow_name, null: false

      json :inputs, null: false
      
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    run('CREATE SEQUENCE figures_ids START 1;')
  end

  down do
    drop_table(:figures)

    run('DROP SEQUENCE figures_ids')
  end
end
