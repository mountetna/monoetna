require 'date'
require 'logger'
require 'rollbar'
require 'tempfile'
require_relative 'helpers'

require 'pry'

class EtnaCommands
  class Help < Etna::Command
    usage 'List this help'

    def execute
      puts 'Commands:'
      EtnaApp.instance.commands.each do |name, cmd|
        puts cmd.usage
      end
    end
  end

  class CreateProject < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'create_project <environment> <filepath>'

    def execute(env, filepath)
      @environ = environment(env)

      create_project_workflow = Etna::Clients::Magma::CreateProjectFromJsonWorkflow.new(
        magma_client: @environ.magma_client,
        janus_client: @environ.janus_client,
        filepath: filepath)
      create_project_workflow.create!
    end
  end

  class ValidateProject < Etna::Command
    include WithLogger
    usage 'validate_project <filepath>'

    def execute(filepath)
      project = Etna::Clients::Magma::JsonProject.new(
        filepath: filepath)
      return puts "Project JSON is well-formatted!" if project.valid?

      puts "Project JSON has #{project.errors.length} errors:"
      project.errors.each do |error|
        puts "  * " + error.gsub("\n", "\n\t")
      end
    end
  end

  class ValidateModel < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'validate_model <environment> <project_name> <model_name> <filepath>'

    def execute(env, project_name, model_name, filepath)
      @environ = environment(env)

      add_model_workflow = Etna::Clients::Magma::AddModelFromJsonWorkflow.new(
        magma_client: @environ.magma_client,
        project_name: project_name,
        model_name: model_name,
        filepath: filepath)
      # If the workflow initialized, then no errors!
      puts "Model JSON for #{model_name} is well-formatted and is valid for project #{project_name}."
    end
  end

  class AddModel < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'add_model <environment> <project_name> <model_name> <filepath>'

    def execute(env, project_name, model_name, filepath)
      @environ = environment(env)

      add_model_workflow = Etna::Clients::Magma::AddModelFromJsonWorkflow.new(
        magma_client: @environ.magma_client,
        project_name: project_name,
        model_name: model_name,
        filepath: filepath)
      add_model_workflow.add!
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Etna instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
    end
  end
end
