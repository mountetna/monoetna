# A vulcan session model.  For now this is just static, but could be made backed by a database, or a browser
# cookie session, or just about anything we like.
class Session
  attr_reader :project_name, :workflow_name, :inputs, :key

  def initialize(project_name, workflow_name, key, inputs = {primary_inputs: {}})
    @project_name = project_name
    @workflow_name = workflow_name
    @key = key
    @inputs = inputs
  end

  def workflow
    @workflow ||= Etna::Workflow.from_yaml_file("#{workflow_name}.cwl")
  end

  def orchestration
    return nil if workflow.nil?
    @orchestration ||= Vulcan::Orchestration.new(workflow, session)
  end

  def define_user_input(step_name, output_name, hash)
    (@inputs[step_name] ||= {})[output_name] = hash
  end

  def include?(step_name)
    @inputs.include?(step_name)
  end

  # Note -- data provided by a user does not, itself, have a script defining how to derive it.
  # We pre-hash it here using a different mechanism than storage cell hashes
  # which depend on a deterministic algorithm and not on the result itself.
  def outputs_hash_for(step_name)
    return nil unless include?(step_name)
    outputs = @inputs[step_name]
    hashes = outputs.keys.sort.map { |output_name| outputs[output_name] }
    Digest::SHA1.hexdigest(hashes.join('&'))
  end

  def output_storage_files(step_name)
    return nil unless include?(step_name)
    outputs = @inputs[step_name]
    ch = outputs_hash_for(step_name)
    return nil if ch.nil?

    outputs.keys.sort.map { |output_name| storage_file_for(step_name, output_name, ch) }
  end

  def storage_file_for(step_name, output_name, ch = nil)
    ch ||= outputs_hash_for(step_name)
    return nil if ch.nil?

    Storage::StorageFile.new(
        project_name: project_name,
        cell_hash: ch,
        data_filename: output_name,
        logical_name: output_name,
    )
  end
end