require 'json'
require 'set'

module Etna
  module Clients
    class Magma
      class ValidatorBase
        attr_reader :errors
        def initialize
          @errors = []
        end

        def valid?
          @errors.length == 0
        end

        def nil_or_empty?(value)
          value.nil? || value.empty?
        end

        def check_key(label, raw, key)
          @errors << "Missing required key for #{label}: \"#{key}\"." if !raw.dig(key)
          @errors << "Invalid empty #{key} for #{label}: \"#{raw[key]}\"." if raw.dig(key) && nil_or_empty?(raw[key])
        end

        def check_in_set(label, raw, key, valid_values)
          @errors << "Invalid #{key} for #{label}: \"#{raw[key]}\".\nShould be one of #{valid_values}." if raw.dig(key) && !valid_values.include?(raw[key])
        end

        def name_regex_with_numbers
          /\A[a-z][a-z0-9]*(_[a-z0-9]+)*\Z/
        end

        def check_valid_name_with_numbers(label, proposed_name)
          @errors << "#{label} name \"#{proposed_name}\" must be snake_case and can only consist of letters, numbers, and \"_\"." unless proposed_name =~ name_regex_with_numbers
        end

        def name_regex_no_numbers
          /\A[a-z]*(_[a-z]+)*\Z/
        end

        def model_exists_in_project?(project_magma_models, model_name)
          !!project_magma_models.model(model_name)
        end
      end

      class ProjectValidator < ValidatorBase
        attr_reader :project
        def initialize(project)
          super()
          @project = project
        end

        def validate
          validate_project_names
          validate_model_links
          validate_models
        end

        def validate_project_names
          check_key('root project', project.raw, 'project_name')
          check_key('root project', project.raw, 'project_name_full')
          name = project.raw['project_name']
          @errors << "Project name #{name} must be snake_case and cannot start with a number or \"pg_\"." unless name =~ name_regex_with_numbers && !name.start_with?('pg_')
        end

        def validate_models
          project.models.all.each do |model|
            validator = ModelValidator.new(model)
            validator.validate
            @errors += validator.errors unless validator.valid?
          end
        end

        def validate_model_links
          # Check all the attributes of type 'link', and make sure
          #   that the link_model_name exists in the project definition.
          model_names = project.models.all.map { |m| m.name }
          project.models.all.map do |model|
            validator = ModelValidator.new(model)
            validator.validate_link_models(model_names)
            @errors += validator.errors unless validator.valid?
          end
        end
      end

      class ModelValidator < ValidatorBase
        attr_reader :model
        def initialize(model)
          super()
          @model = model
          @valid_parent_link_types = Etna::Clients::Magma::ParentLinkType.entries.sort # sort for prettier presentation
        end

        def name
          model.name
        end

        def raw
          model.raw
        end

        def parent_model_name
          model.raw['parent_model_name']
        end

        def is_project?
          model.name == 'project'
        end

        def validate
          validate_add_model_data
          validate_attributes
        end

        def validate_add_model_data
          @errors << "Model name #{name} must be snake_case and can only consist of letters and \"_\"." unless model.name =~ name_regex_no_numbers

          if !is_project?
            check_key("model #{name}", raw, 'parent_model_name')
            check_key("model #{name}", raw, 'parent_link_type')

            check_in_set("model #{name}", raw, 'parent_link_type', @valid_parent_link_types)
          end
          if model.raw['parent_link_type'] != Etna::Clients::Magma::ParentLinkType::TABLE
            check_key("model #{name}", raw, 'identifier')
          end
        end

        def validate_attributes
          model.template.attributes.attribute_keys.each do |attribute_name|
            attribute = model.template.attributes.attribute(attribute_name)
            attribute_validator = AttributeValidator.new(attribute)
            attribute_validator.validate_key(attribute_name)
            attribute_validator.validate_add_attribute_data
            @errors += attribute_validator.errors unless attribute_validator.valid?
          end
        end

        def validate_link_models(model_names)
          @errors << "Parent model \"#{parent_model_name}\" for #{name} does not exist in project.\nCurrent models are #{model_names}." if !is_project? && !model_names.include?(parent_model_name)

          link_attributes.each do |attribute|
            attribute_validator = AttributeValidator.new(attribute)
            attribute_validator.validate_link_models(model_names)
            @errors += attribute_validator.errors unless attribute_validator.valid?

            @errors << "Model \"#{name}\" already belongs to parent model \"#{parent_model_name}\". Remove attribute \"#{attribute.attribute_name}\"." if attribute.link_model_name == parent_model_name
          end
        end

        def validate_existing_models(project_magma_models)
          @errors << "Model #{name} already exists in project!" if model_exists_in_project?(project_magma_models, name)
        end

        def link_attributes
          model.template.attributes.all.select do |attribute|
            attribute.attribute_type == Etna::Clients::Magma::AttributeType::LINK
          end
        end
      end

      class AttributeValidator < ValidatorBase
        attr_reader :attribute
        def initialize(attribute)
          super()
          @attribute = attribute
          # NOTE: for input simplicity, I've removed some of the
          #   link-related types, since we'll try to calculate
          #   those while parsing the JSON structure.
          @valid_attribute_types = Etna::Clients::Magma::AttributeType.entries.reject { |a|
            a == Etna::Clients::Magma::AttributeType::CHILD ||
            a == Etna::Clients::Magma::AttributeType::COLLECTION ||
            a == Etna::Clients::Magma::AttributeType::IDENTIFIER ||
            a == Etna::Clients::Magma::AttributeType::PARENT
          }.sort # sort for prettier presentation

          @valid_validation_types = Etna::Clients::Magma::AttributeValidationType.entries.sort  # sort for prettier presentation
        end

        def validate_basic_attribute_data
          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'attribute_type')

          # The following two could be calculated or left blank?
          # But it would be nice in the final UI to have them be more informative, so
          #   we'll enforce here.
          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'display_name')
          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'desc')

          check_in_set("attribute #{attribute.attribute_name}", attribute.raw, 'attribute_type', @valid_attribute_types)
        end

        def validate_attribute_validation
          return unless attribute.validation
          check_key("attribute #{attribute.attribute_name}, validation", attribute.validation, 'type')
          check_key("attribute #{attribute.attribute_name}, validation", attribute.validation, 'value')
          check_in_set("attribute #{attribute.attribute_name}, validation", attribute.validation, 'type', @valid_validation_types)
        end

        def validate_key(attribute_key)
          @errors << "Attribute key \"#{attribute_key}\" must match attribute_name \"#{attribute.attribute_name}\"." unless attribute_key == attribute.attribute_name
        end

        def is_link_attribute?
          attribute.attribute_type == Etna::Clients::Magma::AttributeType::LINK
        end

        def is_identifier_attribute?
          attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
        end

        def validate_link_models(model_names)
          return unless is_link_attribute?

          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'link_model_name')

          # Check that the linked model exists.
          @errors << "Linked model, \"#{attribute.link_model_name}\", on attribute #{attribute.attribute_name} does not exist!\nCurrent models are #{model_names}." unless model_names.include?(attribute.link_model_name)
        end

        def validate_add_attribute_data
          check_valid_name_with_numbers('Attribute', attribute.attribute_name)
          validate_basic_attribute_data unless is_identifier_attribute?
          validate_attribute_validation

          @errors << "Attribute name #{attribute.attribute_name} should match the link_model_name, \"#{attribute.link_model_name}\"." if is_link_attribute? && attribute.link_model_name != attribute.attribute_name
        end
      end

      class AttributeActionValidatorBase < ValidatorBase
        attr_reader :action, :project_models
        def initialize(action, project_models)
          super()
          @action = action
          @project_models = project_models
        end

        def action_to_attribute(action)
          action_json = JSON.parse(action.to_json)
          # Magma Model uses type and description for Actions,
          #   but Attribute uses attribute_type and desc.
          action_json['attribute_type'] = action_json.delete('type')
          action_json['desc'] = action_json.delete('description')
          Etna::Clients::Magma::Attribute.new(action_json)
        end

        def validate
          raise "Subclasses must implement this method."
        end

        def exists_in_magma_model?(magma_model_name, attribute_name)
          !!project_models.model(magma_model_name).template.attributes.attribute(attribute_name)
        end

        def validate_model_exists(magma_model_name)
          @errors << "Model \"#{magma_model_name}\" does not exist in project." unless model_exists_in_project?(project_models, magma_model_name)
        end

        def check_already_exists_in_model(magma_model_name, attribute_name)
          return unless model_exists_in_project?(project_models, magma_model_name)
          @errors << "Attribute \"#{attribute_name}\" already exists in model #{magma_model_name}." if exists_in_magma_model?(magma_model_name, attribute_name)
        end

        def check_does_not_exist_in_model(magma_model_name, attribute_name)
          return unless model_exists_in_project?(project_models, magma_model_name)
          @errors << "Attribute \"#{attribute_name}\" does not exist in model #{magma_model_name}." unless exists_in_magma_model?(magma_model_name, attribute_name)
        end
      end

      class AddAttributeActionValidator < AttributeActionValidatorBase
        def initialize(action, project_models)
          super
          @attribute = action_to_attribute(action)
        end

        def validate
          validate_attribute_data
          validate_model_exists(action.model_name)
        end

        def validate_attribute_data
          validator = AttributeValidator.new(@attribute)
          validator.validate_basic_attribute_data
          validator.validate_attribute_validation
          @errors += validator.errors unless validator.valid?

          check_already_exists_in_model(action.model_name, action.attribute_name)
        end
      end

      class AddLinkActionValidator < AttributeActionValidatorBase
        def source
          action.links.first
        end

        def dest
          action.links.last
        end

        def validate
          validate_links
          validate_both_models_exist
          validate_link_data
        end

        def validate_links
          check_key("action #{action}", action, :links)
          @errors << "Must include two link entries, each with \"model_name\", \"attribute_name\", and \"type\"." unless action.links.length == 2
          check_key("link #{source}", source, :model_name)
          check_key("link #{source}", source, :attribute_name)
          check_key("link #{source}", source, :type)
          check_key("link #{dest}", dest, :model_name)
          check_key("link #{dest}", dest, :attribute_name)
          check_key("link #{dest}", dest, :type)
        end

        def validate_both_models_exist
          validate_model_exists(source[:model_name])
          validate_model_exists(dest[:model_name])
        end

        def validate_link_data
          # Make sure the attribute names don't already exist in the models,
          #   and that the types are valid.
          link_types = Set.new([source[:type], dest[:type]])
          expected_link_types = Set.new([
            Etna::Clients::Magma::AttributeType::LINK,
            Etna::Clients::Magma::AttributeType::COLLECTION
          ])
          if link_types != expected_link_types
            @errors << "You must have one \"link\" and one \"collection\" type in the links."
          else
            check_already_exists_in_model(source[:model_name], source[:attribute_name])
            check_already_exists_in_model(dest[:model_name], dest[:attribute_name])

            @errors << "Links #{source} and #{dest} must point to each other." unless source[:model_name] == dest[:attribute_name] && source[:attribute_name] == dest[:model_name]
          end
        end
      end

      class RenameAttributeActionValidator < AttributeActionValidatorBase
        def validate
          validate_action
          validate_model_exists(action.model_name)
          validate_proposed_name
        end

        def validate_action
          check_key("action #{action}", action, :model_name)
          check_key("action #{action}", action, :attribute_name)
          check_key("action #{action}", action, :new_attribute_name)
        end

        def validate_proposed_name
          check_does_not_exist_in_model(action.model_name, action.attribute_name)
          check_valid_name_with_numbers('New attribute', action.new_attribute_name)
          check_already_exists_in_model(action.model_name, action.new_attribute_name)
        end
      end

      class UpdateAttributeActionValidator < AttributeActionValidatorBase
        def initialize(action, project_models)
          super
          @attribute = action_to_attribute(action)
        end

        def validate
          validate_attribute_data
          validate_model_exists(action.model_name)
        end

        def validate_attribute_data
          check_does_not_exist_in_model(action.model_name, action.attribute_name)

          validator = AttributeValidator.new(@attribute)
          validator.validate_attribute_validation
          @errors += validator.errors unless validator.valid?
        end
      end

      class AttributeActionsValidator < ValidatorBase
        attr_reader :actions, :project_models
        def initialize(actions, project_models)
          super()
          @actions = actions
          @project_models = project_models
        end

        def project_model_names
          project_models.all.map(&:template).map(&:name)
        end

        def validate
          validate_actions
        end

        def validate_actions
          actions.each do |action|
            clazz = Object.const_get(clazz_name(action.action_name))
            validator = clazz.new(action, project_models)
            validator.validate

            @errors += validator.errors unless validator.valid?
          end
        end

        def camelize(action_name)
          action_name.split('_').map(&:capitalize).join('')
        end

        def clazz_name(action_name)
          "Etna::Clients::Magma::#{camelize(action_name)}ActionValidator"
        end
      end
    end
  end
end
