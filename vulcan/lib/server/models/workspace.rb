class Vulcan
  class Workspace < Sequel::Model
    many_to_one :workflow
    one_to_many :runs
  end

  def to_hash
    {
      id: id,
      workflow_id: workflow_id,
      repo_branch: repo_branch,
      user_email: user_email,
      workspace_dir: workspace_dir,
      created_at: created_at,
      updated_at: updated_at
    }
  end
end
