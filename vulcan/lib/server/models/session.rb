require 'securerandom'

# A vulcan session model.  For now this is just static, but could be made backed by a database, or a browser
# cookie session, or just about anything we like.
class Session < Etna::Cwl
  attr_reader :project_name, :workflow_name, :inputs, :key

  def initialize(attributes)
    @attributes = attributes
    @project_name = attributes['project_name']
    @workflow_name = attributes['workflow_name']
    @workflow_name += ".cwl" unless @workflow_name =~ /\.cwl$/
    @key = attributes['key']
    @inputs = attributes['inputs']
  end

  def self.new_session_for(project_name, workflow_name, key, inputs = {})
    self.new({
        'project_name' => project_name,
        'workflow_name' => workflow_name,
        'key' => key,
        'inputs' => inputs,
    })
  end

  FIELD_LOADERS = {
      project_name: Etna::Cwl::PrimitiveLoader::STRING,
      workflow_name: Etna::Cwl::PrimitiveLoader::STRING,
      key: Etna::Cwl::PrimitiveLoader::STRING.optional.map { |v| v.nil? || v.empty? ? SecureRandom.uuid.hex.to_s : v },
      inputs: Etna::Cwl::StrictMapLoader.new(Etna::Cwl::AnyLoader::ANY.map { |v| {json_payload: v}}, Etna::Cwl::SourceLoader.new).optional.map { |v| v || {} },
  }

  def self.from_json(json)
    json = json.map { |k, v| [ k.to_s, v ] }.to_h
    loader.load(json)
  end

  def self.from_figure(figure)
    self.from_json(figure.to_hash)
  end

  def as_json
    super.tap do |result|
      result['inputs'] = result['inputs'].map do |k, v|
        [Etna::Cwl.source_as_string(k), v[:json_payload]]
      end.to_h
    end
  end

  def run!(storage:, token:)
    orchestration.load_json_inputs!(storage)
    orchestration.scheduler.schedule_more!(orchestration: orchestration, token: token, storage: storage)
  end

  def status(storage:, build_target_cache: {})
    {
      session: self.as_json,
      status: all_steps_status(storage: storage, build_target_cache: build_target_cache),
      outputs: primary_output_status(storage: storage, build_target_cache: build_target_cache)
    }
  end

  def primary_output_status(storage:, build_target_cache:)
    orchestration.build_target_for(:primary_outputs, build_target_cache).tap do |bt|
      return orchestration.scheduler.status(storage: storage, build_target: bt, step: nil).update({
        downloads: bt.is_built?(storage) ? bt.downloads(storage) : nil
      })
    end
  end

  def all_steps_status(storage:, build_target_cache:)
    orchestration.serialized_step_path.map do |paths|
      paths.map do |step_name|
        step = workflow.find_step(step_name)

        next if step.nil?

        step.status_json(
          storage: storage,
          orchestration: orchestration,
          build_target_cache: build_target_cache
        )
      end.compact
    end
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
