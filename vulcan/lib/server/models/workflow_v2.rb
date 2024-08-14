class Vulcan
  class WorkflowV2 < Sequel::Model(:workflows)
    one_to_many :workspaces

    def to_hash
        {
          id: id,
          project: project,
          name: name,
          author: author,
          branch: branch,
          repo_remote_url: repo_remote_url,
          repo_path: repo_path,
          repo_tag: repo_tag,
          created_at: created_at,
          updated_at: updated_at
        }
    end

    def repo_name
      File.basename(repo_path)
    end

  end
end
