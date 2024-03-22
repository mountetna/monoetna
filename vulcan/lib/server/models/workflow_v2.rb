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
          repo_tag: repo_tag,
          config: config,
          created_at: created_at,
          updated_at: updated_at
        }
    end

    def jobs
      config.keys
    end

  end
end
