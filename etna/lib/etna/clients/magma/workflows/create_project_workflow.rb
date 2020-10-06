# Creates a project from a JSON file.
# This workflow:
# 1) Creates the project in Janus.
# 2) Adds the user as a project administrator to Janus.
# 3) Refreshes the user's token with the new privileges.
# 4) Creates the project in Magma.
# 5) Iterates through the models and attributes and creates those in Magma.

require 'json'
require 'ostruct'
require_relative './model_synchronization_workflow'
require_relative '../../janus/models'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class CreateProjectWorkflow < Struct.new(:magma_client, :janus_client, :project_file, keyword_init: true)
        attr_reader :project

        def initialize(**params)
          super({}.update(params))
          @project = Project.new(project_file)
          raise "Project JSON has errors", project.errors unless project.valid?
        end

        def create_janus_project
          janus_client.add_project(Etna::Clients::Janus::AddProjectRequest.new(
            project_name: project.name,
            project_name_full: project.name_full
          ))
        end

        def create!
          create_project_models
          update_model_attributes
        end

        def user
          janus_client.whoami(Etna::Clients::Janus::WhoAmIRequest.new).user
        end
      end

      class JsonProject
        attr_reader :errors
        def initialize(filepath:)
          @raw = JSON.parse(File.read(filepath), symbolize_names: true)
          @errors = []
          validate
        end

        def name
          @raw[:project_name]
        end

        def name_full
          @raw[:project_name_full]
        end

        def models
          @models ||= @raw[:models].map { |model_name, model_def|
            JsonModel.new(model_name, model_def) }
        end

        def validate
          validate_project_names
          validate_models
        end

        def valid?
          errors.length == 0
        end

        def validate_project_names
          errors << "Missing required key for root project, \"project_name\"." if !@raw.key?(:project_name)
          errors << "Missing required key for root project, \"project_name_full\"." if !@raw.key?(:project_name_full)
        end

        def validate_models
          models
        end
      end

      class JsonModel
        attr_reader :name
        def initialize(model_name, raw)
          @name = model_name
          @raw = raw
          validate
        end

        def validate
          validate_add_model_data
          validate_attributes
        end

        def validate_add_model_data
          raise "Missing required key for model #{name}, \"identifier\"." if !@raw.key?(:identifier)
          if name != :project
            raise "Missing required key for model #{name}, \"parent_model_name\"." if !@raw.key?(:parent_model_name)
            raise "Missing required key for model #{name}, \"parent_link_type\"." if !@raw.key?(:parent_link_type)
          end
        end

        def validate_attributes
          # Check that attribute types are supported
          # Check that any validations are supported
          # Check that parent_link_type and link_types are supported
          # Make sure that all parent models are defined
        end
      end
    end
  end
end
