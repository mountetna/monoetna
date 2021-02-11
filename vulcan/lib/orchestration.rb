require 'open3'

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

    def run_until_done!(storage)
      until (runnables = next_runnable_build_targets(storage)).empty?
        run!(storage: storage, build_target: runnables.first)
      end
    end


    def command(storage:, input_files:, output_files:)
      [
          "docker",
          "run",
          "-i",
          "--rm",
          "--mount",
          "-v",
          "/var/run/docker.sock:/var/run/docker.sock:ro",
          "-v",
          "archimedes-exec:/archimedes-exec",
          Vulcan.instance.config(:archimedes_run_image),
          "poetry",
          "run",
          "archimedes-run",
          "--isolator=docker",
          "--image=" + Vulcan.instance.config(:archimedes_image),
      ] + output_files.map do |sf|
        "--output=#{sf.to_archimedes_storage_file(storage)}"
      end + input_files.map do |sf|
        "--input=#{sf.to_archimedes_storage_file(storage)}"
      end
    end

    def run_script!(storage:, script:, input_files:, output_files:)
      status = nil
      output_str = nil

      Open3.popen2(command(storage: storage, input_files: input_files, output_files: output_files)) do |input, output, wait_thr|
        input.write(script)
        output_str = output.read
        status = wait_thr.value.exitstatus
      end

      if status != 0
        raise "Failure to run archimedes-run command: #{output}"
      end

      JSON.parse(output_str)
    end

    def run!(storage:, build_target:)
      storage.with_build_transaction(build_target) do |output_files|
        if build_target.script.is_a?(Hash)
          # This is a copy cell.  Find inputs, link them.
          output_files.each do |of|
            input_file = script[of.logical_name]
            unless input_file.nil?
              ::File.link(input_file.data_path(storage), of.data_path(storage))
            end
          end
          { 'status' => 'done' }
        else
          run_script!(
              storage: storage, input_files: build_target.input_files,
              output_files: output_files, script: build_target.script
          )
        end
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

    def next_runnable_build_targets(storage)
      build_targets_for_paths.map do |path_build_targets|
        path_build_targets.select { |bt| bt.should_build?(storage) }.last
      end.select { |v| !v.nil? }
    end

    def build_targets_for_paths
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

    def take_input_file_for_material_reference(material_reference, input_name)
      Storage::MaterialSource.new(
          project_name: session.project_name, session_key: session.key,
          material_reference: material_reference).take_as_input(input_name)
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
          script = {}
          workflow.inputs.each do |input|
            output_filenames << input.id
            input_files << take_input_file_for_material_reference(
                material_reference_for_user_input([:primary_inputs, input.id], input), input.id)
            script[input.id] = input_files.last
          end
        elsif step_name == :primary_outputs
          workflow.outputs.each do |output|
            output_filenames << output.id

            source_step_name, source_output_name = output.outputSource
            input_files << build_target_for(source_step_name, cache).take_as_input(source_output_name, output.id)
          end
        elsif (step = workflow.find_step(step_name))
          if step.ui_query_name
            script = {}
          elsif step.script_name
            script = step.script
          end

          step.out.each do |step_out|
            output_filenames << step_out.id
            if session.include?(step_name, step_out.id) && step.ui_query_name
              input_files << take_input_file_for_material_reference(session.material_reference_for([step_name, step_out.id]), step_out.id)
              script[step_out.id] = input_files.last
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