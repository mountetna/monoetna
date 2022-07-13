Sequel.migration do
  up do
    # Tests validate these properties
    create_table(:revision_footprints) do
      primary_key :id
      # Not foreign keys -- these need to be resistant to data removal.
      integer :bucket_id, null: false
      # Misnomer -- will use -1 for 'null'
      integer :parent_id, null: false
      # Misnomer -- will use -1 for 'null', due to the way that
      integer :file_id, null: false
      # Misnomer -- will use -1 for 'null', due to the way that
      integer :folder_id, null: false
      DateTime :updated_at, null: false
    end

    alter_table(:revision_footprints) do
      # This is also incidentally an index to make queries easy.
      add_unique_constraint [:bucket_id, :parent_id, :folder_id, :file_id], unique_constraint_name: "unique_revision_footprints"
    end

    run Metis.triggers_sql
  end
end
