Sequel.migration do
  change do
    create_table(:backups) do
      primary_key :id
      String :md5_hash, null: false, unique: true
      String :description, null: false, unique: true
      String :archive_id, null: false
    end
    alter_table(:files) do
      add_foreign_key :backup_id, :backups
    end
    # copy the existing file names into the backup table

    self[:files].exclude(archive_id: nil).all.uniq{|f| f[:file_hash]}.each do |file|
      backup_id = self[:backups].insert(
        md5_hash: file[:file_hash],
        archive_id: file[:archive_id],
        description: "md5:#{file[:file_hash]} #{file[:file_name]}"
      )
      self[:files].where(file_hash: file[:file_hash]).update(backup_id: backup_id)
    end

    alter_table(:files) do
      drop_column :archive_id
    end
  end

  down do
    alter_table(:files) do
      drop_column :backup_id
    end
    drop_table(:backups)
  end
end
