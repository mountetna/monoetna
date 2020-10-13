# Adds a model to a project, from a JSON file.
# This workflow:
# 1) Validates the model JSON file.
# 2) Validates that the model parent and any link attributes exist in Magma.
# 3) Adds the model to the Magma project.
# 4) Adds any attributes to Magma.

require 'json'
require 'ostruct'
require_relative './model_synchronization_workflow'
require_relative './json_models'

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
          @model = Etna::Clients::Magma::JsonModel.new(
            model_name,
            JSON.parse(File.read(filepath)))

          @errors = model.errors.dup
          validate

          raise "Model JSON has errors:\n  * #{format_errors}" unless valid?
        end

        def format_errors
          @errors.map { |e| e.gsub("\n", "\n\t") }.join("\n  * ")
        end

        def valid?
          @errors.length == 0
        end

        def validate
          validate_model_against_project
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
          @errors << "Model #{model.name} already exists in project #{project_name}!" if project_model_names.include?(model.name)

          @errors << "Parent model #{model.parent_model_name} does not exist in project #{project_name}.\nCurrent models are #{project_model_names}." unless project_model_names.include?(model.parent_model_name)

          model.link_attributes do |attribute|
            @errors << "Linked model #{attribute.link_model_name} does not exist in project #{project_name}.\nCurrent models are #{project_model_names}." unless project_model_names.include?(attribute.link_model_name)
          end
        end

        def create_magma_model
          # Technically ensure_magma_tree should create models as needed,
          #   but things seem to go more smoothly if we manually
          #   create the models before we try to ensure the tree.
          magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project_name,
            actions: [Etna::Clients::Magma::AddModelAction.new(
              model_name: model.name,
              parent_model_name: model.parent_model_name,
              parent_link_type: model.parent_link_type,
              identifier: model.identifier
            )]))
        end

        def magma_models
          # Convert our JSON model into an instance of
          #   Magma Model + Attributes
          models = Etna::Clients::Magma::Models.new
          model_builder = models.build_model(model.name)
          model.to_magma_model(model_builder)

          # We also have to include "models" for all linked models, and add in
          #   the reciprocal link attributes.
          model.link_attributes do |attribute|
            model_builder = models.build_model(attribute.link_model_name)
            linked_model = Etna::Clients::Magma::JsonModel.from_name(attribute.link_model_name)
            linked_model.add_reciprocal_link_attribute(model_builder, model)
            linked_model.to_magma_model(model_builder)
          end

          models
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
          magma_models.model(model.name).template.attributes.all.select do |attribute|
            attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
          end.each do |attribute|
            magma_client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
              project_name: project_name,
              actions: [Etna::Clients::Magma::UpdateAttributeAction.new(
                model_name: model.name,
                attribute_name: attribute.name,
                display_name: attribute.display_name,
                description: attribute.desc,
                validation: attribute.validation
              )]))
          end
        end

        def add!
          create_magma_model
          ensure_magma_tree
          update_magma_attributes
        end
      end
    end
  end
end
