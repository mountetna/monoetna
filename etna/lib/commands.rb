require 'date'
require 'logger'
require 'rollbar'
require 'tempfile'
require_relative 'helpers'

module Etna
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

    class AttributeActionsTemplate < Etna::Command
      include WithLogger
      usage 'attribute_actions_template'

      def execute()
        spec = Gem::Specification.find_by_name("etna")
        gem_root = spec.gem_dir
        FileUtils.cp(
          "#{gem_root}/lib/etna/templates/attribute_actions_template.json",
          'attribute_actions_template.json')
        puts "A sample attribute actions JSON template has been provided at `attribute_actions_template.json`."
      end
    end

    class ValidateModel < Etna::Command
      include WithEtnaClientsByEnvironment
      include WithLogger
      usage 'validate_model <environment> <project_name> <model_name> <filepath>'

      def execute(env, project_name, model_name, filepath)
        @environ = environment(env)

        # Use the workflow to validate instead of the ModelValidator directly,
        #   because we also need to do the simple project-level validations, i.e.
        #   does the model exist in the project already?
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

    class ValidateAttributeActions < Etna::Command
      include WithEtnaClientsByEnvironment
      include WithLogger
      usage 'validate_attribute_actions <environment> <project_name> <filepath>'

      def execute(env, project_name, filepath)
        @environ = environment(env)

        # Use the workflow to validate instead of the AttributeActionsValidator directly,
        #   because we also need to check the actions against the project, i.e.
        #   do the right models exist, do the attributes already exist, etc.
        attribute_actions_workflow = Etna::Clients::Magma::AttributeActionsFromJsonWorkflow.new(
          magma_client: @environ.magma_client,
          project_name: project_name,
          filepath: filepath)
        # If the workflow initialized, then no errors!
        puts "Attribute Actions JSON is well-formatted and is valid for project #{project_name}."
      end
    end

    class AttributeActions < Etna::Command
      include WithEtnaClientsByEnvironment
      include WithLogger
      usage 'attribute_actions <environment> <project_name> <filepath>'

      def execute(env, project_name, filepath)
        @environ = environment(env)

        # Use the workflow to validate instead of the AttributeActionsValidator directly,
        #   because we also need to check the actions against the project, i.e.
        #   do the right models exist, do the attributes already exist, etc.
        attribute_actions_workflow = Etna::Clients::Magma::AttributeActionsFromJsonWorkflow.new(
          magma_client: @environ.magma_client,
          project_name: project_name,
          filepath: filepath)
        attribute_actions_workflow.run!
      end
    end

    class UpdateAttributesFromCsv < Etna::Command
      include WithEtnaClientsByEnvironment
      include WithLogger
      usage 'update_attributes_from_csv <environment> <project_name> <model_name> <filepath>'

      def magma_crud
        @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
            magma_client: @environ.magma_client,
            project_name: @project_name)
      end

      def execute(env, project_name, model_name, filepath)
        @environ = environment(env)
        @project_name = project_name

        update_attributes_workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflowSingleModel.new(
            magma_crud: magma_crud,
            project_name: project_name,
            model_name: model_name,
            filepath: filepath)
        update_attributes_workflow.update_attributes
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
end
