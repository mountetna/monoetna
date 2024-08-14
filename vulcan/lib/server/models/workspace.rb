class Vulcan
  class Workspace < Sequel::Model
    many_to_one :workflow_v2, :key => :workflow_id
    one_to_many :runs
  end


  def to_hash
    {
      workspace_id: id,
      workspace_name: name,
      workflow_id: workflow_id,
      workflow_name: workflow.name,
      user_email: user_email,
      workspace_path: path,
      created_at: created_at,
      updated_at: updated_at
    }
  end

  def to_hash_with_vulcan_config
    # Use the basic hash and add vulcan_config
    to_hash.merge(vulcan_config: workflow.vulcan_config)
  end

end
