Sequel.migration do
  change do
    alter_table(:workspaces) do
      add_column :git_request_by, String, null: false, default: "tag_or_sha"
      add_column :git_request, String, null: false, default: ""
    end
  end
end