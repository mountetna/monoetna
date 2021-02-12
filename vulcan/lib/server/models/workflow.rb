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
    end

    class Step < Cwl
      SCRIPT_REGEX = /^scripts\/(.*)$/
      UI_QUERIES_REGEX = /^ui-queries\/(.*)$/

      def lookup_operation_script
        return nil if script_name.nil?

        # TODO: Run some type checking that the operation's inputs / outputs match the step.
        script_path = File.join(Vulcan.instance.config(:workflows_folder), "scripts/#{script_name}.py")
        if ::File.exists?(script_path)
          run.load_attributes_from_script(script_path)
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