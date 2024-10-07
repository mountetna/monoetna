class Vulcan
  class WorkflowV2 < Sequel::Model(:workflows)
    one_to_many :workspaces

    def to_hash
        {
          id: id,
          projects: projects,
          name: name,
          branch: branch,
          repo_remote_url: repo_remote_url,
          created_at: created_at,
          updated_at: updated_at
        }
    end
  end
end
