# Extends the base Etna::Workflow class with extra logic useful for
module Etna
  class Workflow
    def self.from_yaml_file(filename)
      attributes = YAML.safe_load(File.join(File.dirname(__FILE__), "../workflows/#{filename}"))
      self.class.loader.load(attributes)
    end

    def find_step(step_name)
      steps.select { |s| s.id == step_name }
    end

    def find_operation_script(step_name)
      step = self.find_step(step_name)
      return nil if step.nil?
      step.lookup_operation_script
    end
  end

  class Step
    def lookup_operation_script
      run = @attributes['run']
      return nil if run.nil?

      if run.is_a?(Operation)
        # TODO: Run some type checking that the operation's inputs / outputs match the step.
        # Also, would be nice if
        script_path = File.join(File.dirname(__FILE__), "../scripts/#{run.id}.py")
        if ::File.exists?(script_path)
          run.load_attributes_from_script(script_path)
          return ::File.read(script_path)
        end
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