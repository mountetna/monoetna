class Vulcan
  class WorkflowV2 < Sequel::Model(:workflows)
    one_to_many :workspaces

    def to_hash
        {
          id: id,
          project: project,
          # TODO: change to name
          workflow_name: workflow_name,
          author: author,
          repo_remote_url: repo_remote_url,
          repo_local_path: repo_local_path,
          created_at: created_at,
          updated_at: updated_at
        }
    end

    def repo_name
      File.basename(repository_url.to_s)
    end

  end
end
