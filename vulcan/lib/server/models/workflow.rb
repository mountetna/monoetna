# Extends the base Etna::Workflow class with extra logic useful for
module Etna
  class Cwl
    class Workflow < Cwl
      def self.from_yaml_file(filename, prefix = Vulcan.instance.config(:workflows_folder))
        attributes = YAML.safe_load(::File.read(File.join(prefix, filename)))
        self.loader.load(attributes)
      end

      def find_step(step_name)
        steps.find { |s| s.id == step_name }
      end

      def find_operation_script(step_name)
        step = self.find_step(step_name)
        return nil if step.nil?
        step.lookup_operation_script
      end

      def as_steps_json
        {
            class: "Workflow",
            cwlVersion: @attributes['cwlVersion'],
            inputs: @attributes['inputs'].map(&:as_steps_inputs_json_pair).to_h,
            outputs: @attributes['outputs'].map(&:as_steps_output_json_pair).to_h,
            steps: Vulcan::Orchestration.unique_paths(self).map do |path|
              path.map { |step_name| self.find_step(step_name)&.as_step_json }.select { |v| v }
            end
        }
      end
    end

    class WorkflowInputParameter < Cwl
      def as_steps_inputs_json_pair
        [self.id, {
            label: @attributes['label'],
            type: @attributes['type'],
            format: @attributes['format'],
            default: @attributes['default'],
        }]
      end
    end

    class WorkflowOutputParameter < Cwl
      def as_steps_output_json_pair
        [self.id, {
            outputSource: output_source_as_string,
            label: @attributes['label'],
            type: @attributes['type'],
            default: @attributes['default'],
            format: @attributes['format'],
        }]
      end

      def output_source_as_string
        os = @attributes['outputSource']
        if os.first.is_a?(Symbol)
          return os.last
        end

        "#{os.first}/#{os.last}"
      end
    end

    # TODO: Add type checking capabilities that unfold the run operation or workflow and checks that the given
    # in and out match sanely to the operation or workflow typings.
    class Step < Cwl
      SCRIPT_REGEX = /^scripts\/(.*)\.cwl$/
      UI_QUERIES_REGEX = /^ui-queries\/(.*)$/

      def as_step_json
        {
            name: @attributes['id'],
            run: @attributes['run'].id,
            in: @attributes['in'].map(&:id),
            out: @attributes['out'].map(&:id),
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
        run = @attributes['run']
        return nil if run.nil?

        if run.is_a?(Operation)
          m = SCRIPT_REGEX.match(run.id)
          m.nil? ? nil : m[1]
        end
      end

      def ui_query_name
        run = @attributes['run']
        return nil if run.nil?

        if run.is_a?(Operation)
          m = UI_QUERIES_REGEX.match(run.id)
          m.nil? ? nil : m[1]
        end
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