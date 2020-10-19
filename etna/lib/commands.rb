require 'date'
require 'logger'
require 'rollbar'
require 'tempfile'
require_relative 'helpers'

class EtnaApp
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
      user_json = JSON.parse(File.read(filepath))
      magma_json = Etna::Clients::Magma::ProjectConverter.convert_project_user_json_to_magma_json(user_json)
      @project = Etna::Clients::Magma::Project.new(magma_json)

      @validator = Etna::Clients::Magma::ProjectValidator.new(@project)
      @validator.validate

      return puts "Project JSON is well-formatted!" if @validator.valid?

      puts "Project JSON has #{@validator.errors.length} errors:"
      @validator.errors.each do |error|
        puts "  * " + error.gsub("\n", "\n\t")
      end
    end
  end

  class ProjectTemplate < Etna::Command
    include WithLogger
    usage 'project_template'

    def execute()
      spec = Gem::Specification.find_by_name("etna")
      gem_root = spec.gem_dir
      FileUtils.cp(
        "#{gem_root}/lib/etna/templates/create_project_template.json",
        'project_template.json')
      puts "A sample project JSON template has been provided at `project_template.json`."
    end
  end

  class ModelTemplate < Etna::Command
    include WithLogger
    usage 'model_template'

    def execute()
      spec = Gem::Specification.find_by_name("etna")
      gem_root = spec.gem_dir
      FileUtils.cp(
        "#{gem_root}/lib/etna/templates/add_model_template.json",
        'model_template.json')
      puts "A sample model JSON template has been provided at `model_template.json`."
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
