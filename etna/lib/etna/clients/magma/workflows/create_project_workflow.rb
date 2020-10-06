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
          raise "Project JSON has errors: ", project.errors unless project.valid?
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

      class JsonBase
        attr_reader :errors
        def initialize
          @errors = []
        end

        def valid?
          errors.length == 0
        end

        def nil_or_empty?(value)
          value.nil? || value.empty?
        end

        def check_key(label, raw, key)
          errors << "Missing required key for #{label}: \"#{key}\"." if !raw.key?(key.to_sym)
          errors << "Invalid #{key} for #{label}: \"#{raw[key.to_sym]}\"." if raw.key?(key.to_sym) && nil_or_empty?(raw[key.to_sym])
        end

        def check_type(label, raw, key, valid_types)
          errors << "Invalid #{key} for #{label}: \"#{raw[key.to_sym]}\"." if raw.key?(key.to_sym) && !valid_types.include?(raw[key.to_sym])
        end
      end

      class JsonProject < JsonBase
        def initialize(filepath:)
          super()
          @raw = JSON.parse(File.read(filepath), symbolize_names: true)
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

          # Make sure that all parent models are defined
        end

        def validate_project_names
          check_key('root project', @raw, :project_name)
          check_key('root project', @raw, :project_name_full)
        end

        def validate_models
          models.each { |model|
            @errors += model.errors unless model.valid?
          }
        end
      end

      class JsonModel < JsonBase
        attr_reader :name
        def initialize(model_name, raw)
          super()
          @name = model_name
          @raw = raw

          @valid_parent_link_types = [
            Etna::Clients::Magma::ParentLinkType::CHILD,
            Etna::Clients::Magma::ParentLinkType::COLLECTION,
            Etna::Clients::Magma::ParentLinkType::TABLE
          ]
          validate
        end

        def attributes
          @attributes ||= @raw.key?(:attributes) ? @raw[:attributes].map { |attribute_name, attribute_def|
            JsonAttribute.new(attribute_name, attribute_def) } : []
        end

        def validate
          validate_add_model_data
          validate_attributes
        end

        def validate_add_model_data
          check_key("model #{name}", @raw, :identifier)
          if name != :project
            check_key("model #{name}", @raw, :parent_model_name)
            check_key("model #{name}", @raw, :parent_link_type)
            check_type("model #{name}", @raw, :parent_link_type, @valid_parent_link_types)
          end
        end

        def validate_attributes
          return unless @raw.key?(:attributes)

          # Check that attribute types are supported
          # Check that any validations are supported
          # Check that link_types are supported
        end
      end

      class JsonAttribute < JsonBase
        attr_reader :name
        def initialize(attribute_name, raw)
          super()
          @name = attribute_name
          @raw = raw
          validate
        end

        def validate
        end
      end
    end
  end
end
