require 'ostruct'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class ModelSynchronizationWorkflow < Struct.new(:target_client, :source_models, :target_project, keyword_init: true)
        def target_models
          @target_models ||= begin
            target_client.retrieve(RetrievalRequest.new(project_name: self.target_project, model_name: 'all')).models
          end
        end

        def execute_updates(*actions)
          update = UpdateModelRequest.new(project_name: self.target_project)
          actions.each { |a| update.add_action(a) }
          @target_models = nil
          target_client.update_model(update)
        end

        def self.from_api_source(source_project:, source_client:, **kwds)
          self.new(
              source_models: source_client.retrieve(RetrievalRequest.new(project_name: source_project, model_name: 'all')).models,
              **kwds
          )
        end

        # Subclass and override when the source <-> target mapping is not 1:1.
        def target_of_source(model_name)
          model_name
        end

        # Subclass and override when the source <-> target attribute mapping is not 1:1.
        def target_attribute_of_source(model_name, attribute_name)
          attribute_name
        end

        # Potentially cyclical, protected against re-entry via the seen cache.
        # Establishes the link attributes in a given model graph.
        def ensure_model_tree(model_name, seen = Set.new)
          return unless (source_model = source_models.model(model_name))
          return if seen.include?(model_name)
          seen.add(model_name)
          ensure_model(model_name)

          attributes = source_model.template.attributes

          attributes.all.each do |attribute|
            # Don't copy or update parent links.  Once a model has been setup with a parent someway.
            unless attribute.attribute_type == AttributeType::PARENT
              if attribute.link_model_name
                ensure_model(attribute.link_model_name)
                ensure_model_link(model_name, attribute.link_model_name, attribute.attribute_name)
              else
                ensure_model_attribute(model_name, attribute.attribute_name)
              end
            end


            # Even if it's a parent node, however, we still want to cascade the tree expansion to all links.
            unless attribute.link_model_name.nil?
              ensure_model_tree(attribute.link_model_name, seen)
            end
          end
        end

        def ensure_model_link(model_name, link_model_name, attribute_name)
          return unless (model = source_models.model(model_name))
          return unless (source_attribute = model.template.attributes.attribute(attribute_name))

          return unless (link_model = source_models.model(link_model_name))
          link_model_attributes = link_model.template.attributes
          reciprocal = link_model_attributes.all.find do |attr|
            attr.link_model_name == model_name
          end

          target_model_name = target_of_source(model_name)
          target_link_model_name = target_of_source(link_model_name)

          target_attributes = target_models.model(target_model_name).template.attributes
          return if target_attributes.attribute_keys.include?(attribute_name)

          add_link = AddLinkAction.new
          add_link.links << AddLinkDefinition.new(model_name: target_model_name, attribute_name: attribute_name, type: source_attribute.attribute_type)
          add_link.links << AddLinkDefinition.new(model_name: target_link_model_name, attribute_name: reciprocal.attribute_name, type: reciprocal.attribute_type)

          execute_updates(add_link)
        end

        def ensure_model_attribute(model_name, attribute_name)
          return unless (model = source_models.model(model_name))
          return unless (source_attribute = model.template.attributes.attribute(attribute_name))

          target_model_name = target_of_source(model_name)
          target_attribute_name = target_attribute_of_source(model_name, attribute_name)

          target_attributes = target_models.model(target_model_name).template.attributes
          return if target_attributes.attribute_keys.include?(attribute_name)

          add_attribute = AddAttributeAction.new(
              model_name: target_model_name,
              attribute_name: target_attribute_name,
          )

          Attribute.copy(source_attribute, add_attribute)

          execute_updates(add_attribute)
        end

        # Non cyclical, non re-entrant due to the requirement that parents cannot form a cycle.
        # This method, and it's partner prepare_parent, should never call into any re-entrant or potentially
        # cyclical method, like ensure_model_tree.
        def ensure_model(model_name)
          return unless (source_model = source_models.model(model_name))

          target_model_name = target_of_source(model_name)
          return if target_models.model_keys.include?(target_model_name)

          template = source_model.template

          add_model_action = AddModelAction.new(
              {
                  model_name: target_model_name,
                  identifier: template.identifier,
              }
          )

          parents = template.attributes.all.select { |a| a.attribute_type == AttributeType::PARENT }
          parent_model_name, parent_link_type = prepare_parent(model_name, template, parents)
          unless parent_model_name.nil?
            add_model_action.parent_model_name = parent_model_name
            add_model_action.parent_link_type = parent_link_type
          end

          execute_updates(add_model_action)
        end

        # Non cyclical, non re-entrant due to the requirement that parents cannot form a cycle.
        # This method, and it's partner and ensure_model, should never call into any re-entrant or potentially
        # cyclical method, like ensure_model_tree.
        def prepare_parent(model_name, template, parents)
          return [nil, nil] if parents.empty?
          return [nil, nil] unless (parent_model = source_models.model(template.parent))
          return [nil, nil] unless (child_attribute = parent_model.template.attributes.attribute(model_name))

          ensure_model(template.parent)

          [
              target_of_source(template.parent),
              child_attribute.attribute_type
          ]
        end
      end

      class ShallowCopyModelWorkflow < ModelSynchronizationWorkflow
        def initialize(model_name:, **kwds)
          super(**kwds)
          @model_name = model_name
        end

        # Aside from just creating models for links, do not cascade the expansion.
        def ensure_model_tree(model_name, *args)
          puts "Checking #{model_name} #{@model_name}"
          if model_name == @model_name
            super(model_name, *args)
          end
        end
      end
    end
  end
end
