Sequel.migration do
  up do
    create_table(:workflow_snapshots) do
      primary_key :id

      foreign_key :figure_id, :figures, null: false, unique: true
      json :authors, null: false
      json :projects
      String :vignette, text: true
      String :display_name
      String :description
      String :cwl_yaml, null: false
      String :query_action
      json :input_query_map
      json :tags
      String :image
      json :scripts, null: false

      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:workflow_snapshots)
  end
end
