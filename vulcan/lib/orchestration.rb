class Vulcan
  # This is the temporary glue for running synchronous workflows from a given session and its inputs.
  # Most of this logic will evolve outside of vulcan depending on where we land on orchestration.
  class Orchestration
    attr_reader :workflow, :session
    def initialize(workflow, session)
      @workflow = workflow
      @session = session
    end

    # Returns a list of lists, describing the dependency between steps, primary inputs, and primary outputs.
    # Each inner list is a unique path in side of the workflow starting a primary_inputs and terminating
    # at either a primary_output or a step that is not used as an input the workflow
    def unique_paths
      @unique_paths ||= begin
        directed_graph = ::DirectedGraph.new

        workflow.steps.each do |step|
          step.in.each do |step_input|
            directed_graph.add_connection(step_input.source.first, step.id)
          end
        end

        workflow.outputs.each do |output|
          directed_graph.add_connection(output.outputSource.first, :primary_outputs)
        end

        directed_graph.paths_from(:primary_inputs)
      end
    end


    def self.find_runnable_steps
      unique_paths.map do |path|
        path.find do |step_name|
          step = workflow.find_step(step_name)
          input_files = input_files_for(step)
          output_files = output_files_for(step)

          input_files.all? { |input_file| ::File.exists?(input_file.data_path) } &&
              output_files.any? { |output_file| ! ::File.exists?(output_file.data_path) }
        end
      end
    end

    def input_files_for(step)
      step.in.map { |step_input| create_input_storage_file(step_input.id, step_input.source) }
    end

    def output_files_for(step)
      step.out.map { |step_output| create_output_storage_file(step.id, step_output.id) }
    end

    def cell_hash(step_name)
      if session.include?(step)
        return session.outputs_hash_for(step)
      end

      step = workflow.find_step(step_name)
      return nil if step.nil?

      input_files = []
      step.in.each do |step_input|
        inner_step, inner_output_name = step_input.source
        ch = cell_hash(inner_step)
        return nil if ch.nil?

        input_files << Storage::StorageFile.new(
            project_name: session.project_name, cell_hash: ch,
            data_filename: inner_output_name, logical_name: step_input.id
        )
      end

      script = workflow.find_operation_script(step.id)
      return nil if script.nil?

      Storage.cell_hash(
          project_name: session.project_name,
          input_files: input_files,
          output_filenames: step.out.map(&:id),
          session_key: session.key,
          script: script,
      )
    end

    def self.create_input_storage_file(input_name, source)
      step_name, output_file = source
      ch = cell_hash(step_name)

      return nil if ch.nil?

      Storage::StorageFile.new(
          project_name: session.project_name, cell_hash: ch,
          data_filename: output_file, logical_name: input_name
      )
    end

    def create_output_storage_file(step_name, output_name)
      ch = cell_hash(step_name)

      return nil if ch.nil?

      Storage::StorageFile.new(
          project_name: session.project_name, cell_hash: ch,
          data_filename: output_name, logical_name: output_name,
      )
    end
  end
end