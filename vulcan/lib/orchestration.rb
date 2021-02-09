class Vulcan
  # This is the temporary glue for running synchronous workflows from a given session and its inputs.
  # Most of this logic will evolve outside of vulcan depending on where we land on orchestration.
  class Orchestration
    attr_reader :workflow, :session

    def initialize(workflow, session)
      @workflow = workflow
      @session = session
    end

    MAX_RUNNABLE = 20

    def run!
      i = 0
      # Because this is synchronous right now, and we want to prevent any bad code causing infinite loop
      # that consumes all system resources, we restrict ourselves to a maximum of 20 iterations to complete work.
      # In a future async sysystem this restriction can go away.
      until (runnable = find_runnable_steps.find { |s| !s.nil? }).nil? or (i += 1) > MAX_RUNNABLE
        step = workflow.find_step(runnable)
        # input_files = input_files_for(step)
        # output_files = output_files_for(step)
      end
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


    def find_build_targets
      build_target_cache = {}

      unique_paths.map do |path|
        path.map do |step_name|
          step = workflow.find_step(step_name)
          build_target_for(step.id, build_target_cache)
        end
      end
    end

    # Combines any session given value for a primary input and the potential default that may be defined for it.
    def material_reference_for_user_input(source, input)
      if session.include?(source)
        session.material_reference_for(source)
      elsif input.default
        {json_payload: JSON.dump(input.default)}
      else
        {unfulfilled: source}
      end
    end

    # Cache is used to save effort in computing common recursive targets within a single calculation session.
    # Do not pass an actual value in for it.
    def build_target_for(step_name, cache = {})
      if cache.include?(step_name)
        return cache[step_name]
      end

      cache[step_name] = begin
        output_filenames = []
        input_files = []
        script = nil

        if step_name == :primary_inputs
          workflow.inputs.each do |input|
            output_filenames << input.id
            input_files << Storage::MaterialSource.new(
                project_name: session.project_name, session_key: session.key,
                material_reference: material_reference_for_user_input([:primary_inputs, input.id], input))
          end
        elsif step_name == :primary_outputs
          workflow.outputs.each do |output|
            output_filenames << output.id

            source_step_name, source_output_name = output.outputSource
            input_files << build_target_for(source_step_name, cache).take_as_input(source_output_name, output.id)
          end
        elsif (step = workflow.find_step(step_name))
          step.out.each do |step_out|
            output_filenames << step_out.id
            if session.include?(step_name, step_out.id)
              input_files << session.material_reference_for([step_name, step_out.id]).as_input(step_out.id)
            end
          end

          step.in.each do |step_in|
            source_step_name, source_output_name = step_in.source
            input_files << build_target_for(source_step_name, cache).take_as_input(source_output_name, step_in.id)
          end
        end

        Storage::BuildTarget.new(
            project_name: session.project_name,
            session_key: session.key,
            input_files: input_files,
            output_filenames: output_filenames,
            script: script,
        )
      end
    end
  end
end