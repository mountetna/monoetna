Sequel.migration do
  change do
    alter_table(:workspaces) do
      add_column :git_ref, String, null: false, default: ""
    end
    
    alter_table(:workspaces) do
      rename_column :git_version, :git_sha
    end
  end
end