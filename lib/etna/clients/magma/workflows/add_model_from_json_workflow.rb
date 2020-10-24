# Adds a model to a project, from a JSON file.
# This workflow:
# 1) Validates the model JSON file.
# 2) Validates that the model parent and any link attributes exist in Magma.
# 3) Adds the model to the Magma project.
# 4) Adds any attributes to Magma.

require 'json'
require 'ostruct'
require_relative './model_synchronization_workflow'
require_relative './json_validators'
require_relative './json_converters'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class AddModelFromJsonWorkflow < Struct.new(:magma_client, :project_name, :model_name, :filepath, keyword_init: true)
        attr_reader :model
        def initialize(**params)
          super({}.update(params))
          user_json = JSON.parse(File.read(filepath))

          magma_json = Etna::Clients::Magma::ModelConverter.convert_model_user_json_to_magma_json(
            model_name, user_json)

          @model = Etna::Clients::Magma::Model.new(magma_json)

          @validator = Etna::Clients::Magma::ModelValidator.new(model)
          @validator.validate

          validate_model_against_project

          raise "Model JSON has errors:\n  * #{format_errors(@validator.errors)}" unless @validator.valid?

          @converter = ModelConverter.new(model)
          @converter.convert!
        end

        def format_errors(errors)
          errors.map { |e| e.gsub("\n", "\n\t") }.join("\n  * ")
        end

        def project_models
          @project_models ||= magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
            project_name: project_name,
            model_name: 'all')).models
        end

        def project_model_names
          project_models.all.map(&:template).map(&:name)
        end

        def validate_model_against_project
          # Make sure that the parent model exists.
          # Make sure any linked models exist.
          @validator.validate_existing_models(project_models)
          @validator.validate_link_models(project_model_names)
        end

        def create_magma_model
          # Technically ensure_magma_tree should create models as needed,
          #   but things seem to go more smoothly if we manually
          #   create the models before we try to ensure the tree.
          magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project_name,
            actions: [Etna::Clients::Magma::AddModelAction.new(
              model_name: model.name,
              parent_model_name: model.raw['parent_model_name'],
              parent_link_type: model.raw['parent_link_type'],
              identifier: model.template.identifier
            )]))
        end

        def magma_models
          all_models = project_models + model

          model.template.attributes.all.select do |attribute|
            attribute.attribute_type == Etna::Clients::Magma::AttributeType::LINK
          end.each do |attribute|
            linked_model = all_models.model(attribute.link_model_name)
            link_converter = Etna::Clients::Magma::ModelConverter.new(linked_model)
            link_converter.add_reciprocal_link_attribute(model)
          end

          all_models
        end

        def ensure_magma_tree
          # Call the model_synchronization_workflow.ensure_model_tree
          #   on each model to set their attributes in the target Magma.
          workflow = Etna::Clients::Magma::ModelSynchronizationWorkflow.new(
            target_project: project_name,
            target_client: magma_client,
            source_models: magma_models
          ).ensure_model_tree(model.name)
        end

        def update_magma_attributes
          # Sometimes an identifier attribute might have
          #   validations or other settings included in the
          #   model definition. So we have to "update" those
          #   instead of adding those (which occurs in
          #   ensure_magma_tree).
          model.template.attributes.all.select do |attribute|
            attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
          end.each do |attribute|
            magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
              project_name: project_name,
              actions: [Etna::Clients::Magma::UpdateAttributeAction.new(
                model_name: model.name,
                attribute_name: attribute.attribute_name,
                display_name: attribute.display_name,
                description: attribute.desc,
                validation: attribute.validation
              )]))
          end
        end

        def add!
          puts "Adding the new model to Magma."
          create_magma_model
          ensure_magma_tree
          update_magma_attributes
          puts "All complete!"
        end
      end
    end
  end
end
