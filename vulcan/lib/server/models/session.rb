# A vulcan session model.  For now this is just static, but could be made backed by a database, or a browser
# cookie session, or just about anything we like.
class Session
  attr_reader :project_name, :workflow_name, :inputs

  def initialize(project_name, workflow_name, inputs)
    @project_name = project_name
    @workflow_name = workflow_name
    @inputs = inputs
  end
end