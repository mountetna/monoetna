Sequel.migration do
  change do
    create_table(:figures) do
      Integer :id, null: false

      String :title

      json :document, null: false
      
      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end

    run('CREATE SEQUENCE figures_ids START 1;')
  end
end
