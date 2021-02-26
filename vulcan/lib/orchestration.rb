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
      load_json_inputs!(storage)
      [].tap do |ran|
        errors = RunErrors.new
        i = 0
        until (runnables = next_runnable_build_targets(storage)).empty? || (i += 1) > MAX_RUNNABLE
          begin
            next_runnable = runnables.find { |r| !errors.include?(r) }
            break if next_runnable.nil?

            run!(storage: storage, build_target: runnables.first)
          rescue => e
            errors << {cell_hash: runnables.first.cell_hash, error: e}
          ensure
            ran << runnables.first
          end
        end
        raise errors unless errors.empty?
      end
    end

    def load_json_inputs!(storage)
      workflow.inputs.each do |input|
        reference = material_reference_for_user_input([:primary_inputs, input.id], input)
        load_json_payload!(storage, reference)
      end

      session.material_references.each do |reference|
        load_json_payload!(storage, reference)
      end
    end

    # Note: do not validate json payloads here -- do so either at the entry point or simply allow cells to fail with
    # bad inputs.  Validating here would start to conflate potential storage or logic errors with simple user input
    # errors.
    def load_json_payload!(storage, reference)
      if (payload = reference[:json_payload])
        ms = material_source(reference)
        unless ::File.exists?(ms.build_output&.data_path(storage))
          storage.with_build_transaction(ms) do |build_files|
            path = build_files.first.data_path(storage)
            ::File.write(path, payload)
          end
        end
      end
    end

    def command(storage:, input_files:, output_files:)
      [
          "docker",
          "run",
          "-i",
          "--rm",
          "-v",
          "/var/run/docker.sock:/var/run/docker.sock:ro",
          "-v",
          "#{Vulcan.instance.config(:archimedes_exec_volume)}:/archimedes-exec",
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

      cmd = command(storage: storage, input_files: input_files, output_files: output_files)
      Open3.popen2(*cmd) do |input, output, wait_thr|
        input.print(script)
        input.close
        output_str = output.read
        status = wait_thr.value.exitstatus
      end

      if status != 0
        Vulcan.instance.logger.warn("archimedes-run failure: #{output_str}")
        raise "Failure to run archimedes-run command!  Check logs."
      end

      JSON.parse(output_str)
    end

    def relative_path(from, to)
      Pathname.new(to).relative_path_from(Pathname.new(::File.dirname(from))).to_s
    end

    def run!(storage:, build_target:)
      storage.with_build_transaction(build_target) do |output_files|
        script = build_target.script
        if script.is_a?(Hash)
          # This is a copy cell.  Find inputs, link them.
          output_files.each do |of|
            output_path = of.data_path(storage)
            input_file = script[of.logical_name]
            unless input_file.nil?
              input_path = input_file.data_path(storage)
              relative_path = relative_path(output_path, input_path)
              # puts("coping for #{of}")
              # puts(output_path, input_path, relative_path)
              ::File.unlink(output_path)
              ::File.symlink(relative_path, output_path)
            end
          end
          {'status' => 'done'}
        elsif script
          result = run_script!(
              storage: storage, input_files: build_target.input_files,
              output_files: output_files, script: script
          )

          if (error = result["error"])
            raise "Python error while executing script: #{error}"
          end
        else
          raise "Could not determine or find backing script"
        end
      end

      raise "Build for #{build_target} failed to produce outputs!" unless build_target.is_built?(storage)
    end

    # Returns a list of lists, describing the dependency between steps, primary inputs, and primary outputs.
    # The first inner list is a serialized path traversing all nodes in a strictly linear fashion.
    # Note that this ordering is NOT Guaranteed to be the execution order; tasks that exist on divergent paths
    # might be run in parallel, and thus may complete in a non deterministic ordering.  But this serialized
    # ordering will, nonetheless, occur in an order such that an strict dependency is always before
    # convergent nodes.
    # After the first inner list, each other inner list is a unique path in side of the workflow starting
    # a primary_inputs and terminating at either a primary_output or a step that is not used as an input the workflow
    # Each path represents a potential divergent -> convergent ordering that has strict execution ordering.  The
    # first node and last node of each of these paths _may_ be shared with other paths as part of convergence.
    def self.unique_paths(workflow)
      directed_graph = ::DirectedGraph.new

      workflow.steps.each do |step|
        step.in.each do |step_input|
          directed_graph.add_connection(step_input.source.first, step.id)
        end
      end

      workflow.outputs.each do |output|
        directed_graph.add_connection(output.outputSource.first, :primary_outputs)
      end

      [ directed_graph.serialized_path_from(:primary_inputs) ] + directed_graph.paths_from(:primary_inputs)
    end

    def unique_paths
      @unique_paths ||= self.class.unique_paths(workflow)
    end

    def next_runnable_build_targets(storage)
      build_targets_for_paths.map do |path_build_targets|
        path_build_targets.select { |bt|
          bt.should_build?(storage)
        }.last
      end.select { |v| !v.nil? }.uniq { |bt| bt.cell_hash }
    end

    def build_targets_for_paths
      build_target_cache = {}

      unique_paths.map do |path|
        path.map do |step_name|
          step = workflow.find_step(step_name)
          build_target_for(step&.id || step_name, build_target_cache)
        end
      end
    end

    def primary_input_material_sources
      workflow.inputs.map do |input|
        material_source(material_reference_for_user_input([:primary_inputs, input.id], input))
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
          script = {}
          workflow.inputs.zip(primary_input_material_sources).each do |input, source|
            output_filenames << input.id
            input_files << source.take_as_input(input.id)
            script[input.id] = input_files.last
          end
        elsif step_name == :primary_outputs
          script = {}
          workflow.outputs.each do |output|
            output_filenames << output.id

            source_step_name, source_output_name = output.outputSource
            input_file = build_target_for(source_step_name, cache).take_as_input(source_output_name, output.id)
            script[output.id] = input_file
            if input_file.nil?
              raise "Could not find output #{source_output_name.inspect} from step #{source_step_name.inspect} for primary output #{output.id.inspect}"
            end
            input_files << input_file
          end
        elsif (step = workflow.find_step(step_name))
          if step.ui_query_name
            script = {}
          elsif step.script_name
            script = step.lookup_operation_script
            raise "Could not find backing script #{step.script_name.inspect} for step #{step.id}" if script.nil?
          else
            raise "Step #{step.id} has invalid run: #{step.run}.  Must be either a ui-queries/ or scripts/ entry." if script.nil?
          end

          step.out.each do |step_out|
            output_filenames << step_out.id
            if step.ui_query_name
              ref = session.material_reference_for([step_name, step_out.id])
              source = material_source(ref)
              input_files << source.take_as_input(step_out.id)
              script[step_out.id] = input_files.last
            end
          end

          step.in.each do |step_in|
            source_step_name, source_output_name = step_in.source
            input_file = build_target_for(source_step_name, cache).take_as_input(source_output_name, step_in.id)
            if input_file.nil?
              raise "Could not find output #{source_output_name.inspect} from source #{source_step_name.inspect} while building input #{step_in.id.inspect} for step #{step.id.inspect}"
            end
            input_files << input_file
          end
        else
          raise "Step #{step_name} has no backing definition."
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

    def material_source(material_reference)
      Storage::MaterialSource.new(
          project_name: session.project_name, session_key: session.key,
          material_reference: material_reference)
    end

    class RunErrors < StandardError
      def initialize
        @errors = {}
      end

      def <<(value)
        @errors[value[:cell_hash]] = value[:error]
      end

      def include?(bt)
        @errors.key?(bt.cell_hash)
      end

      def message_for_build_target(bt)
        @errors[bt.cell_hash].message
      end

      def message
        @errors
      end

      def empty?
        @errors.keys.empty?
      end
    end
  end
end