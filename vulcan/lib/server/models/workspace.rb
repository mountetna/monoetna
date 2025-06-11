class Vulcan
  class Workspace < Sequel::Model
    many_to_one :workflow_v2, :key => :workflow_id
    one_to_many :configs

    def to_hash
      {
        workspace_id: id,
        name: name,
        workflow_id: workflow_id,
        workflow_name: workflow_v2.name,
        user_email: user_email,
        tags: tags,
        workspace_path: path,
        git_ref: git_ref,
        git_sha: git_sha,
        created_at: created_at,
        updated_at: updated_at
      }
    end
  end
end
