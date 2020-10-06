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
          validate_model_links
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

        def validate_model_links
          # Check all the attributes of type 'link', and make sure
          #   that the link_model_name exists in the project definition.
          model_names = models.map { |m| m.name }

          models.each { |model|
            model.attributes.each { |attribute|
              if attribute.type == Etna::Clients::Magma::AttributeType::LINK
                check_key("model #{model.name}, attribute #{name}", attribute.raw, :link_model_name)

                # Check that the linked model exists.
                errors << "Model \"#{model.name}\" already belongs to parent model \"#{model.parent_model_name}\". Remove attribute \"#{attribute.name}\"." if attribute.link_model_name == model.parent_model_name

                errors << "Linked model, \"#{attribute.link_model_name}\", on attribute #{attribute.name} of model #{model.name} does not exist!" if !model_names.include?(attribute.link_model_name)
              end
            }
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

        def parent_model_name
          @raw[:parent_model_name]&.to_sym
        end

        def identifier
          @raw[:identifier]&.to_sym
        end

        def attributes
          @attributes ||= @raw.key?(:attributes) ? @raw[:attributes].map { |attribute_name, attribute_def|
            JsonAttribute.new(self, attribute_name, attribute_def) } : []
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
          attributes.each { |attribute|
            @errors += attribute.errors unless attribute.valid?
          }
        end
      end

      class JsonAttribute < JsonBase
        attr_reader :name, :model, :raw
        def initialize(model, attribute_name, raw)
          super()
          @model = model
          @name = attribute_name
          @raw = raw

          # NOTE: for input simplicity, I've removed some of the
          #   link-related types, since we'll try to calculate
          #   those while parsing the JSON structure.
          @valid_attribute_types = [
            Etna::Clients::Magma::AttributeType::STRING,
            Etna::Clients::Magma::AttributeType::DATE_TIME,
            Etna::Clients::Magma::AttributeType::BOOLEAN,
            Etna::Clients::Magma::AttributeType::FILE,
            Etna::Clients::Magma::AttributeType::FLOAT,
            Etna::Clients::Magma::AttributeType::IMAGE,
            Etna::Clients::Magma::AttributeType::INTEGER,
            Etna::Clients::Magma::AttributeType::LINK,
            Etna::Clients::Magma::AttributeType::MATCH,
            Etna::Clients::Magma::AttributeType::MATRIX,
            Etna::Clients::Magma::AttributeType::TABLE,
          ]

          @valid_validation_types = [
            Etna::Clients::Magma::AttributeValidationType::REGEXP,
            Etna::Clients::Magma::AttributeValidationType::ARRAY,
            Etna::Clients::Magma::AttributeValidationType::RANGE
          ]

          # We will calculate this?
          # @valid_link_types = [
          #   Etna::Clients::Magma::AttributeType::COLLECTION,
          #   Etna::Clients::Magma::AttributeType::LINK
          # ]

          validate
        end

        def type
          @raw[:attribute_type]
        end

        def link_model_name
          @raw[:link_model_name]&.to_sym
        end

        def validate
          validate_add_attribute_data
        end

        def validate_add_attribute_data
          check_key("model #{model.name}, attribute #{name}", @raw, :attribute_name)

          if model.identifier != name
            check_key("model #{model.name}, attribute #{name}", @raw, :attribute_type)

            # The following two could be calculated or left blank?
            # But it would be nice in the final UI to have them be more informative, so
            #   we'll enforce here.
            check_key("model #{model.name}, attribute #{name}", @raw, :display_name)
            check_key("model #{model.name}, attribute #{name}", @raw, :desc)

            check_type("model #{model.name}, attribute #{name}", @raw, :attribute_type, @valid_attribute_types)
          end

          if @raw.key?(:validation)
            check_key("model #{model.name}, attribute #{name}, validation", @raw[:validation], :type)
            check_key("model #{model.name}, attribute #{name}, validation", @raw[:validation], :value)
            check_type("model #{model.name}, attribute #{name}, validation", @raw[:validation], :type, @valid_validation_types)
          end
        end
      end
    end
  end
end
