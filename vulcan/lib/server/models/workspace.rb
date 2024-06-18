class Vulcan
  class Workspace < Sequel::Model
    many_to_one :workflow_v2, :key => :workflow_id
    one_to_many :runs
  end

  def to_hash
    {
      id: id,
      workflow_id: workflow_id,
      workflow_name: workflow.name,
      user_email: user_email,
      path: path,
      created_at: created_at,
      updated_at: updated_at
    }
  end
end
