Sequel.migration do
  change do
    create_table(:search_modalities) do
      primary_key :id
      String :name, null: false
      String :model, null: false
      Integer :position, null: false, default: 0
      DateTime :created_at
      DateTime :updated_at

      index [:position, :id]
      index Sequel.function(:lower, :name), unique: true, name: :search_modalities_lower_name_index
    end

    create_table(:search_fields) do
      primary_key :id
      String :label, null: false
      column :query, :json, null: false
      String :kind, null: false
      foreign_key :modality_id, :search_modalities, null: true, on_delete: :restrict
      Integer :position, null: false, default: 0
      DateTime :created_at
      DateTime :updated_at

      check(kind: %w[count filter])

      index [:position, :id]
      index :modality_id
    end
  end
end
