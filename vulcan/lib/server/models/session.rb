# A vulcan session model.  For now this is just static, but could be made backed by a database, or a browser
# cookie session, or just about anything we like.
class Session
  attr_reader :project_name, :workflow_name, :inputs, :key

  def initialize(project_name, workflow_name, key, inputs = {})
    @project_name = project_name
    @workflow_name = workflow_name
    @key = key
    @inputs = inputs
  end

  def as_json
    {
        project_name: project_name,
        key: key,
        inputs: inputs,
        workflow_name: workflow_name,
    }
  end

  def self.from_json(json)
    self.new(json['project_name'], json['workflow_name'], json['key'], json['inputs'])
  end

  def workflow
    @workflow ||= Etna::Cwl::Workflow.from_yaml_file(workflow_name)
  end

  def orchestration
    return nil if workflow.nil?
    @orchestration ||= Vulcan::Orchestration.new(workflow, self)
  end

  def define_user_input(source, json_obj)
    @inputs[source] = {json_payload: JSON.dump(json_obj)}
  end

  def include?(source)
    @inputs.include?(source)
  end

  def material_references
    @inputs.values
  end

  def material_reference_for(source)
    if include?(source)
      @inputs[source]
    else
      {unfulfilled: source}
    end
  end
end