Sequel.migration do
  up do
    run <<-SQL
      CREATE SEQUENCE IF NOT EXISTS revision_updates_version_seq;
    SQL

    create_table(:revision_updates) do
      primary_key :id
      # Not foreign keys -- these need to be resistant to data removal.
      integer :bucket_id, null: false
      integer :folder_id, null: true
      integer :file_id, null: true
      integer :version, null: false
      boolean :moved, null: false
    end

    run Metis.triggers_sql
  end
end
