# Extends the base Etna::Workflow class with extra logic useful for
module Etna
  class Cwl
    def self.source_as_string(os)
      if os.first.is_a?(Symbol)
        return os.last
      end

      "#{os.first}/#{os.last}"
    end

    class WorkflowNotFound < Exception
    end

    class Workflow < Cwl
      def self.from_yaml_file(filename, prefix = Vulcan.instance.config(:workflows_folder))
        full_path = File.join(prefix, filename)
        raise WorkflowNotFound.new("File not found: #{filename}") unless ::File.exists?(full_path)
        attributes = YAML.safe_load(::File.read(full_path))
        self.loader.load(attributes)
      end

      def self.metadata(filename, prefix = Vulcan.instance.config(:workflows_folder))
        metadata_file = ::File.join(
          prefix,
          filename.sub("cwl", "metadata.json")
        )

        return {} unless ::File.exists?(metadata_file)

        JSON.parse(::File.read(metadata_file), symbolize_names: true)
      end

      def find_step(step_name)
        steps.find { |s| s.id == step_name }
      end

      def find_operation_script(step_name)
        step = self.find_step(step_name)
        return nil if step.nil?
        step.lookup_operation_script
      end

      def as_steps_json(name)
        {
          class: "Workflow",
          cwlVersion: @attributes["cwlVersion"],
          name: name,
          inputs: @attributes["inputs"].map(&:as_steps_inputs_json_pair).to_h,
          outputs: @attributes["outputs"].map(&:as_steps_output_json_pair).to_h,
          steps: Vulcan::Orchestration.serialize_step_path(self).map do |path|
            path.map { |step_name| self.find_step(step_name)&.as_step_json }.select { |v| v }
          end,
        }
      end

      def step_key(step)
        step_name, var_name = step.source
        if step_name == :primary_inputs
          return :"primary_inputs_#{var_name}"
        end
        return step_name
      end

      def step_graph
        graph = ::DirectedGraph.new.tap do |directed_graph|
          inputs.each do |input|
            directed_graph.add_connection(:root, step_key(input))
          end

          steps.each do |step|
            if step.in.empty?
              directed_graph.add_connection(:root, step.id)
            end

            step.in.each do |step_input|
              directed_graph.add_connection(step_key(step_input), step.id)
            end
          end

          outputs.each do |output|
            directed_graph.add_connection(output.outputSource.first, :primary_outputs)
          end
        end

        graph
      end

      def source_as_string(source)
        Etna::Cwl.source_as_string(source)
      end
    end

    class WorkflowInputParameter < Cwl
      def source
        [:primary_inputs, self.id]
      end

      def as_steps_inputs_json_pair
        [self.id, {
          label: @attributes["label"],
          type: @attributes["type"],
          format: @attributes["format"],
          default: @attributes["default"],
          doc: @attributes["doc"],
        }]
      end
    end

    class WorkflowOutputParameter < Cwl
      def as_steps_output_json_pair
        [self.id, {
          outputSource: output_source_as_string,
          label: @attributes["label"],
          type: @attributes["type"],
          default: @attributes["default"],
          format: @attributes["format"],
        }]
      end

      def output_source_as_string
        Etna::Cwl.source_as_string(@attributes["outputSource"])
      end
    end

    # TODO: Add type checking capabilities that unfold the run operation or workflow and checks that the given
    # in and out match sanely to the operation or workflow typings.
    class Step < Cwl
      SCRIPT_REGEX = /^scripts\/(.*)\.cwl$/
      UI_QUERIES_REGEX = /^ui-queries\/(.*)$/
      UI_OUTPUTS_REGEX = /^ui-outputs\/(.*)$/

      def as_step_json
        {
          name: @attributes["id"],
          run: @attributes["run"].id,
          in: @attributes["in"].map { |i| input_as_json(i) },
          out: @attributes["out"].map(&:id),
          label: @attributes["label"],
          doc: @attributes["doc"],
        }
      end

      def status_json(storage:, build_target_cache:, orchestration:)
        bt = orchestration.build_target_for(id, build_target_cache)

        orchestration.scheduler.status(storage: storage, build_target: bt, step: self).update({
          name: id,
          hash: bt.cell_hash,
          downloads: (!ui_output? && bt.is_built?(storage)) ? bt.downloads(storage) : nil,
        })
      end

      def input_as_json(input)
        {
          id: input.id,
          source: Etna::Cwl.source_as_string(input.source),
        }
      end

      def lookup_operation_script
        return nil if script_name.nil?

        script_path = File.join(Vulcan.instance.config(:workflows_folder), "scripts/#{script_name}.py")
        if ::File.exists?(script_path)
          ::File.read(script_path)
        end
      end

      def script_name
        run = @attributes["run"]
        return nil if run.nil?

        if run.is_a?(Operation)
          m = SCRIPT_REGEX.match(run.id)
          m.nil? ? nil : m[1]
        end
      end

      def ui_query_name
        run = @attributes["run"]
        return nil if run.nil?

        if run.is_a?(Operation)
          m = UI_QUERIES_REGEX.match(run.id)
          m.nil? ? nil : m[1]
        end
      end

      def ui_output_name
        run = @attributes["run"]
        return nil if run.nil?

        if run.is_a?(Operation)
          m = UI_OUTPUTS_REGEX.match(run.id)
          m.nil? ? nil : m[1]
        end
      end

      def ui_behavior?
        ui_query_name || ui_output?
      end

      def ui_output?
        !!ui_output_name
      end
    end

    class Operation
      # Would be interesting to use python type inference to determine input / output types from a script file.
      def load_attributes_from_script(script_path)
        # TODO
      end
    end
  end
end
