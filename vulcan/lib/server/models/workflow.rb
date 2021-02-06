# A vulcan workflow model.  For now this is just static, but could be made backed by a database, or just about anything we like.
# TODO: We really need a
class Workflow
  attr_reader :workflow_name, :cwl, :description, :display_name

  def initialize(workflow_name, cwl, description: "", display_name: "")
    @workflow_name = workflow_name
    @cwl = cwl
    @description = description
    @display_name = display_name
  end

  # Returns a list of lists, describing the dependency between steps, primary inputs, and primary outputs.
  # Each inner list is a unique path in side of the workflow starting a primary_inputs and terminating
  # at either a primary_output or a step that is not used as an input the workflow
  def unique_paths
    @ordered_steps ||= [].tap do |result|
      directed_graph = ::DirectedGraph.new

      cwl.steps.each do |step|
        step.in.each do |k, ref|
          directed_graph.add_connection(ref.first, step.step_name)
        end
      end

      primary_outputs.keys.sort.each do |output_name|
        ref = primary_outputs[output_name]
        directed_graph.add_connection(ref, :primary_outputs)
      end

      directed_graph.paths_from(:primary_inputs)
    end
  end

  def find_step(step_name)
    cwl.steps.select { |s| s.step_name == step_name }
  end

  # Resolves a string of the forms "a-primary-identifier" or "step_name/output" into
  #   [:primary_inputs, "a-primary-identifier"] or
  #   ["step_name", "output"] respectively
  def self.resolve_reference(reference)
    parts = reference.split('/', max=2)
    if parts.length == 1
      [:primary_inputs, reference]
    end
      parts
    nil
  end

  # Collects primary inputs of the format { name: [step_name, output_filename] }
  def primary_outputs
    @primary_outputs ||= {}.tap do |result|
      outputs = @raw['outputs'] || {}
      outputs.each do |k, input_def|
        unless (src = input_def['outputSource']).nil? || (type = input_def['type']).nil?
          if type == 'File'
            result[k] = Workflow.resolve_reference(src)
          end
        end
      end
    end
  end

  # Collects primary inputs from the cwl definition in the form { name: default_value }
  def primary_inputs
    @primary_inputs ||= {}.tap do |result|
      inputs = @raw['inputs'] || {}

      # TODO: More robust handling for all CWL types.
      inputs.each do |k, input_def|
        unless (default = input_def['default']).nil? || (type = input_def['type']).nil?
          case type
          when 'int'
            value = default.to_i
          when 'boolean'
            case default
            when TrueClass
              value = default
            when FalseClass
              value = default
            when "true"
              value = true
            else
              value = false
            end
          else
            next
          end

          result[k] = value
        end
      end
    end
  end

  class Script
    attr_accessor :script_name, :contents
    def initialize(script_name:, contents:)
      @script_name = script_name
      @contents = contents
    end
  end
end