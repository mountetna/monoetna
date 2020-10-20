# Updates, renames, or adds attributes to a model, from a JSON file.
# This workflow:
# 1) Validates the attributes JSON file.
# 2) Validates that the model(s) and any link attribute models exist in Magma.
# 3) Executes each of the updates. This will have to be direct Magma client calls,
#     since the model_synch workflow only does Adds and does not currently
#     handle other attribute actions.

require 'json'
require 'ostruct'
require_relative './json_models'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class AttributeActionsJsonWorkflow < Struct.new(:magma_client, :project_name, :filepath, keyword_init: true)

        def initialize(**params)
          super({}.update(params))

          actions_json = JSON.parse(File.read(filepath))

          converter = Etna::Clients::Magma::AttributeActionsConverter.new(actions_json)
          actions = converter.convert

          validator = Etna::Clients::Magma::AttributeActionsValidator.new(
            actions,
            project_models)
          validator.validate

          raise "Attributes JSON has errors:\n  * #{format_errors(validator.errors)}" unless validator.valid?
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

        def validate_attribute_names_unique

        end

        def check_required_key(action, key, options)
          @errors << "Action requires \"#{key}\". Options are: #{options}." unless action.key?(key)
        end

        def validate_models_exist
          # Make sure that the parent model exists.
          # Make sure any linked models exist.
          available_actions = ['add_attribute', 'rename_attribute', 'update_attribute']
          required_keys = ['action_name', 'model_name', 'attribute_name']
          @raw.each do |attribute_action|
            check_required_key(attribute_action, 'action_name', available_actions)
            check_required_key(attribute_action, 'model_name', project_model_names)
            check_required_key(attribute_action, 'attribute_name', project_models.model(attribute_action['model_name']).template.attributes.all.map(&.attribute_name)) unless attribute_action['action_name'] == 'add_attribute'
          end
          @errors << "Model does not #{model.name} already exists in project #{project_name}!" if project_model_names.include?(model.name)

          model.validate_link_models(project_model_names)
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
            linked_model = Etna::Clients::Magma::JsonModel.linking_stub_from_name(attribute.link_model_name)
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

        def run!
          puts "Executing the attribute actions against Magma."
          create_magma_model
          ensure_magma_tree
          update_magma_attributes
          puts "All complete!"
        end
      end
    end
  end
end
