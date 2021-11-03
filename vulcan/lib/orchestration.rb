require 'open3'
require 'benchmark'
require_relative './asynchronous_scheduler'

class Vulcan
  # This is the temporary glue for running synchronous workflows from a given session and its inputs.
  # Most of this logic will evolve outside of vulcan depending on where we land on orchestration.
  class Orchestration
    attr_reader :workflow, :session

    # Global threaded scheduler
    @@scheduler = AsynchronousScheduler.new

    def initialize(workflow, session)
      @workflow = workflow
      @session = session
    end

    MAX_RUNNABLE = 20

    def scheduler
      @@scheduler
    end

    # Synchronous scheduler that runs all the processes inline.
    def run_until_done!(storage, token = nil)
      load_json_inputs!(storage)
      [].tap do |ran|
        errors = RunErrors.new
        i = 0
        until (runnables = next_runnable_build_targets(storage)).empty? || (i += 1) > MAX_RUNNABLE
          next_runnable = runnables.find { |(step_name, r)| !errors.include?(r) }&.last

          begin
            break if next_runnable.nil?

            run!(storage: storage, build_target: next_runnable, token: token)
          rescue => e
            errors << { cell_hash: next_runnable.cell_hash, error: e }
          ensure
            ran << next_runnable
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

      # Ensure that the primary inputs, at least, have been loaded so that the status endpoint can report
      # meaningful steps running.
      workflow.inputs.each do |input|
        pi = build_target_for(:primary_inputs, var_name: input.id)
        if pi.should_build?(storage)
          run!(storage: storage, build_target: pi, token: nil)
        end
      end
    end

    # Note: do not validate json payloads here -- do so either at the entry point or simply allow cells to fail with
    # bad inputs.  Validating here would start to conflate potential storage or logic errors with simple user input
    # errors.
    def load_json_payload!(storage, reference)
      if (reference.key?(:json_payload))
        # We want to write `false` files as well,
        #   so do not use payload in the if conditional.
        payload = reference[:json_payload]

        # But we don't want to write `nil` values...
        return if nil == payload

        ms = material_source(reference)
        unless ::File.exists?(ms.build_output&.data_path(storage))
          storage.with_run_cell_build(run_cell: storage.run_cell_for(ms)) do |build_dir, build_files|
            path = build_files.first.data_path(storage)
            # If the payload is a Hash, we need to save it as a standard JSON
            #   string, not a Ruby Hash with symbolized keys.
            ::File.write(path, payload.is_a?(Hash) ? payload.to_json : payload)
            storage.install_build_output(buildable: ms, build_dir: build_dir)
          end
        end
      end
    end

    def command(storage:, input_files:, output_files:, token:, ch:)
      [
        "docker",
        "run",
        "--rm",
      ] + docker_run_args(storage: storage, input_files: input_files, output_files: output_files, token: token, ch: ch)
    end

    def docker_run_args(storage:, input_files:, output_files:, token:, ch:)
      [
        "-i",
        "-v",
        "/var/run/docker.sock:/var/run/docker.sock:ro",
        "-v",
        "#{Vulcan.instance.config(:archimedes_exec_volume)}:/archimedes-exec",
        Vulcan.instance.config(:archimedes_run_image),
        "poetry",
        "run",
        "archimedes-run",
        "--isolator=docker",
        "-e",
        "MAGMA_HOST=#{Vulcan.instance.config(:magma)&.dig(:host)}",
        "-e",
        "TOKEN=#{token}",
        "-e",
        "PROJECT_NAME=#{session.project_name}",
        "--image=" + Vulcan.instance.config(:archimedes_run_image),
      ] + output_files.map do |sf|
        "--output=#{sf.to_archimedes_storage_file(storage)}"
      end + input_files.map do |sf|
        "--input=#{sf.to_archimedes_storage_file(storage)}"
      end
    end

    def run_script!(storage:, script:, input_files:, output_files:, token:, ch:)
      status = nil
      output_str = nil

      cmd = command(
        storage: storage,
        input_files: input_files,
        output_files: output_files,
        token: token,
        ch: ch,
      )
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

    def run_as_copy_cell(storage:, build_target:, script:)
      if script.is_a?(Hash)
        storage.with_run_cell_build(run_cell: storage.run_cell_for(build_target)) do |build_dir, output_files|
          # This is a copy cell.  Find inputs, link them.
          output_files.each do |of|
            output_path = of.data_path(storage)
            input_file = script[of.logical_name]

            unless input_file.nil?
              input_path = input_file.data_path(storage)
              relative_path = relative_path(output_path, input_path)
              ::File.unlink(output_path)
              ::File.symlink(relative_path, output_path)
            end
          end

          storage.install_build_output(buildable: build_target, build_dir: build_dir)
        end

        return true
      end

      false
    end

    # This version keeps a lock on the storage for the entire duration of the script.
    def run_script_synchronously(storage:, build_target:, script:, token:)
      storage.with_run_cell_build(run_cell: storage.run_cell_for(build_target)) do |build_dir, output_files|
        result = run_script!(
          storage: storage,
          input_files: build_target.input_files,
          output_files: output_files,
          script: script,
          token: token,
          ch: build_target.cell_hash,
        )

        if (error = result["error"])
          raise "Python error while executing script: #{error}"
        end

        storage.install_build_output(
          buildable: build_target,
          build_dir: build_dir
        )
      end
    end

    # Synchronous runner implementation that expects the transaction to complete.
    def run!(storage:, build_target:, token:)
      Yabeda.vulcan.job_runtime.measure({script_hash: Storage.hash_json_obj(build_target.script)}, Benchmark.realtime do
        script = build_target.script
        if run_as_copy_cell(storage: storage, build_target: build_target, script: script)
          { 'status' => 'done' }
        elsif script
          run_script_synchronously(
            storage: storage,
            build_target: build_target,
            script: script,
            token: token
          )
        else
          raise "Could not determine or find backing script"
        end

        raise "Build for #{build_target} failed to produce outputs!" unless build_target.is_built?(storage)
      end)
    end

    def unique_paths
      @unique_paths = workflow.step_graph.paths_from(:root, false)
    end

    def self.serialize_step_path(workflow)
      [workflow.step_graph.serialized_path_from(:root, false)]
    end

    def serialized_step_path
      @serialized_step_path ||= self.class.serialize_step_path(workflow)
    end

    def next_runnable_build_targets(storage)
      build_targets_for_paths.map do |path_build_targets|
        path_build_targets.select { |step_name,bt|
          bt.should_build?(storage)
        }.to_a.last
      end.compact.uniq { |step_name,bt| bt.cell_hash }
    end

    def build_targets_for_paths
      build_target_cache = {}

      unique_paths.map do |path|
        path.map do |step_name|
          [ step_name, build_target_for(step_name, build_target_cache) ]
        end.to_h
      end
    end

    def primary_input_material_sources
      workflow.inputs.map do |input|
        material_source(material_reference_for_user_input([:primary_inputs, input.id], input))
      end
    end

    # Cache is used to save effort in computing common recursive targets within a single calculation session.
    # Do not pass an actual value in for it.
    def build_target_for(step_name, cache = {}, var_name: nil)
      cache_name = "#{step_name}#{var_name}"
      if cache.include?(cache_name)
        return cache[cache_name]
      end

      cache[cache_name] = begin
        output_filenames = []
        input_files = []
        script = nil

        if step_name == :primary_inputs
          script = {}
          workflow.inputs.zip(primary_input_material_sources).each do |input, source|
            next if var_name && var_name != input.id

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
          if step.ui_behavior?
            script = {}
          elsif step.script_name
            script = step.lookup_operation_script
            raise "Could not find backing script #{step.script_name.inspect} for step #{step.id}" if script.nil?
          else
            raise "Step #{step.id} has invalid run: #{step.run}.  Must be either a ui-queries/, ui-outputs/, or scripts/ entry." if script.nil?
          end

          step.out.each do |step_out|
            output_filenames << step_out.id
            if step.ui_behavior?
              ref = session.material_reference_for([step_name, step_out.id])
              source = material_source(ref)
              input_files << source.take_as_input(step_out.id)
              script[step_out.id] = input_files.last
            end
          end

          step.in.each do |step_in|
            source_step_name, source_output_name = step_in.source
            input_file = build_target_for(source_step_name, cache, var_name: source_output_name).take_as_input(source_output_name, step_in.id)
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
      elsif nil != input.default
        { json_payload: JSON.dump(input.default) }
      else
        { unfulfilled: source }
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
