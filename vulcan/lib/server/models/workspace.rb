class Vulcan
  class Workspace < Sequel::Model
    many_to_one :workflow
    one_to_many :runs
    def_column_alias :hash_value, :hash
  end

  def to_hash
    {
      id: id,
      workflow_id: workflow_id,
      repo_branch: repo_branch,
      user_email: user_email,
      # TODO: probably can just change this to dir
      workspace_dir: workspace_dir,
      hash: hash_value,
      created_at: created_at,
      updated_at: updated_at
    }
  end
end
