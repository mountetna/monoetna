Sequel.migration do
  up do
    create_table(:workflow_snapshots) do
      primary_key :id

      foreign_key :figure_id, :figures, null: false, unique: true
      json :authors, null: false
      String :name, null: false
      json :projects
      String :vignette, text: true
      String :display_name
      String :description
      String :cwl_version, null: false, default: "v1.1"
      json :inputs
      json :outputs
      json :steps, null: false
      json :other_metadata
      json :tags
      String :image

      DateTime :created_at, null: false
      DateTime :updated_at, null: false
    end
  end

  down do
    drop_table(:workflow_snapshots)
  end
end
