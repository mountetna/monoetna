# A vulcan workflow model.  For now this is just static, but could be made backed by a database, or just about anything we like.
module Etna
  class Workflow
    def self.from_yaml_file(filename)
      attributes = YAML.safe_load(File.join(File.dirname(__FILE__), "../workflows/#{filename}"))
      self.class.loader.load(attributes)
    end

    # Returns a list of lists, describing the dependency between steps, primary inputs, and primary outputs.
    # Each inner list is a unique path in side of the workflow starting a primary_inputs and terminating
    # at either a primary_output or a step that is not used as an input the workflow
    def unique_paths
      @ordered_steps ||= begin
        directed_graph = ::DirectedGraph.new

        steps.each do |step|
          step.in.each do |step_input|
            directed_graph.add_connection(step_input.source.first, step.id)
          end
        end

        outputs.each do |output|
          directed_graph.add_connection(output.outputSource.first, :primary_outputs)
        end

        directed_graph.paths_from(:primary_inputs)
      end
    end

    def find_step(step_name)
      steps.select { |s| s.id == step_name }
    end
  end
end