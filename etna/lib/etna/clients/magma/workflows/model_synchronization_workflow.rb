require 'ostruct'

module Etna
  module Clients
    class Magma
      # Note!  These workflows are not perfectly atomic, nor perfectly synchronized due to nature of the backend.
      # These primitives are best effort locally synchronized, but cannot defend the backend or simultaneous
      # system updates.
      class ModelSynchronizationWorkflow < Struct.new(
          :target_client, :source_models, :target_project, :update_block,
          :model_name, :plan_only, :validate, :use_versions, :renames,
          keyword_init: true)
        def target_models
          @target_models ||= begin
            target_client.retrieve(RetrievalRequest.new(project_name: self.target_project, model_name: 'all')).models
          end
        end

        def planned_actions
          @planned_actions ||= []
        end

        def queue_update(action)
          if plan_only
            plan_update(action)
          else
            execute_update(action)
          end
        end

        def execute_planned!
          @planned_actions.each do |action|
            execute_update(action)
          end
        end

        def execute_update(action)
          update = UpdateModelRequest.new(project_name: self.target_project)
          update_block.call(action) if update_block
          update.add_action(action)
          @target_models = nil
          target_client.update_model(update)
        end

        def copy_link_into_target(link, reciprocal)
          attr = target_models.build_model(link.model_name).build_template.build_attributes.build_attribute(link.attribute_name)
          attr.attribute_type = link.type
          attr.name = attr.attribute_name = link.attribute_name
          attr.link_model_name = reciprocal.model_name
        end

        # Applies the given action to the source models and 'plans' its execution.
        def plan_update(action)
          case action
          when UpdateAttributeAction
            attribute_update = target_models.build_model(action.model_name).build_template.build_attributes.build_attribute(action.attribute_name)
            Attribute.copy(action, attribute_update, attributes: Attribute::EDITABLE_ATTRIBUTE_ATTRIBUTES)
          when AddAttributeAction
            new_attribute = target_models.build_model(action.model_name).build_template.build_attributes.build_attribute(action.attribute_name)
            Attribute.copy(action, new_attribute)
          when AddLinkAction
            first_link = action.links[0]
            second_link = action.links[1]
            copy_link_into_target(first_link, second_link)
            copy_link_into_target(second_link, first_link)
          when AddModelAction
            template = target_models.build_model(action.model_name).build_template
            template.name = action.model_name
            template.identifier = action.identifier
            template.parent = action.parent_model_name

            child_link = AddLinkDefinition.new(type: AttributeType::PARENT, model_name: template.name, attribute_name: template.parent)
            parent_link = AddLinkDefinition.new(type: action.parent_link_type, model_name: template.parent, attribute_name: template.name)

            copy_link_into_target(child_link, parent_link)
            copy_link_into_target(parent_link, child_link)

            ['created_at', 'updated_at'].each do |time_attr_name|
              template.build_attributes.build_attribute(time_attr_name).tap do |attr|
                attr.attribute_type = Etna::Clients::Magma::AttributeType::DATE_TIME
                attr.attribute_name = time_attr_name
              end
            end

            if action.parent_link_type != Etna::Clients::Magma::AttributeType::TABLE
              template.build_attributes.build_attribute(template.identifier).tap do |attr|
                attr.attribute_name = template.identifier
                attr.attribute_type = Etna::Clients::Magma::AttributeType::IDENTIFIER
              end
            end
          when RenameAttributeAction
            attributes = target_models.model(action.model_name).template.attributes
            attributes.raw[action.new_attribute_name] = attributes.raw.delete(action.attribute_name)
          else
            raise "Unexpected plan_action #{action}"
          end

          planned_actions << action
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
          return if target_attributes.attribute_keys.include?(target_link_model_name)

          add_link = AddLinkAction.new
          add_link.links << AddLinkDefinition.new(model_name: target_model_name, attribute_name: attribute_name, type: source_attribute.attribute_type)
          add_link.links << AddLinkDefinition.new(model_name: target_link_model_name, attribute_name: reciprocal.attribute_name, type: reciprocal.attribute_type)

          queue_update(add_link)
        end

        def ensure_model_attribute(model_name, attribute_name)
          return unless (model = source_models.model(model_name))
          return unless (source_attribute = model.template.attributes.attribute(attribute_name))

          target_model_name = target_of_source(model_name)
          target_attribute, target_attribute_name = ensure_model_attribute_target_rename(model_name, attribute_name)

          if target_attribute.nil?
            add_attribute = AddAttributeAction.new(
                model_name: target_model_name,
                attribute_name: target_attribute_name,
            )

            Attribute.copy(source_attribute, add_attribute)
            queue_update(add_attribute)
          else
            # If there are is no diff, don't produce an action.
            target_attribute = Attribute.new(target_attribute.raw)
            target_attribute.set_field_defaults!

            source_attribute = Attribute.new(source_attribute.raw)
            source_attribute.set_field_defaults!

            if !source_attribute.is_edited?(target_attribute)
              return
            end

            update_attribute = UpdateAttributeAction.new(
                model_name: target_model_name,
                attribute_name: target_attribute_name,
            )

            Attribute.copy(source_attribute, update_attribute, attributes: Attribute::EDITABLE_ATTRIBUTE_ATTRIBUTES)
            queue_update(update_attribute)
          end
        end

        # Returns a tuple of the target's attribute, post rename if necessary, if it exists, and the name of the target attribute
        # cases here:
        # 1.  There is no rename for the attribute
        #  a. There target attribute already exists -> [target_attribute, attribute_name]
        #  b. The target attribute does not exist -> [nil, attribute_name]
        # 2.  There is an expected rename from the source
        #  a. The target has neither the renamed attribute or the original attribute -> [nil, new_attribute_name]
        #  b. The target has the renamed attribute already -> [renamed_attribute, new_attribute_name]
        #  c. The target has the source attribute, which is not yet renamed. -> [renamed_attribute, new_attribute_name]
        def ensure_model_attribute_target_rename(model_name, attribute_name)
          target_model_name = target_of_source(model_name)
          target_attribute_name = target_attribute_of_source(model_name, attribute_name)
          return nil unless (target_model = target_models.model(target_model_name))

          target_original_attribute = target_model.template.attributes.attribute(target_attribute_name)

          if renames && (attribute_renames = renames[model_name]) && (new_name = attribute_renames[attribute_name])
            new_name = target_attribute_of_source(model_name, new_name)

            unless target_model.template.attributes.attribute_keys.include?(new_name)
              if target_original_attribute
                rename = RenameAttributeAction.new(model_name: target_model_name, attribute_name: target_attribute_name, new_attribute_name: new_name)
                queue_update(rename)
              else
                return [nil, new_name]
              end
            end

            return [target_model.template.attributes.attribute(new_name), new_name]
          end

          [target_original_attribute, target_attribute_name]
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

          queue_update(add_model_action)
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

        def self.models_affected_by(action)
          case action
          when Etna::Clients::Magma::RenameAttributeAction
            [action.model_name]
          when Etna::Clients::Magma::UpdateAttributeAction
            [action.model_name]
          when Etna::Clients::Magma::AddAttributeAction
            [action.model_name]
          when Etna::Clients::Magma::AddModelAction
            [action.model_name]
          when Etna::Clients::Magma::AddLinkAction
            action.links.map(&:model_name)
          else
            []
          end
        end
      end
    end
  end
end
