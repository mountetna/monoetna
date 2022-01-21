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

        def check_key_empty(label, raw, key)
          @errors << "Invalid key for #{label}: \"#{key}\"." if raw.dig(key)
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

        def validate!(err_message)
          @errors = []
          self.validate
          raise "#{err_message}\n#{format_errors}" unless valid?
        end

        def format_errors
          errors.map { |e| e.gsub("\n", "\n\t") }.join("\n  * ")
        end
      end

      class ProjectValidator < ValidatorBase
        attr_reader :create_args

        def initialize(**create_args)
          super()
          @create_args = create_args
        end

        def validate
          validate_project_names
        end

        def validate_project_names
          check_key('root project', create_args, :project_name)
          check_key('root project', create_args, :project_name_full)
          name = create_args[:project_name]
          @errors << "Project name #{name} must be snake_case and cannot start with a number or \"pg_\"." unless name =~ name_regex_with_numbers && !name.start_with?('pg_')
        end
      end

      class RenamesValidator < ValidatorBase
        attr_reader :models, :renames
        def initialize(models = Models.new, renames = {})
          @models = models
          @renames = renames
          super()
        end

        def validate
          renames.each do |model_name, attribute_renames|
            attribute_renames.each do |old_name, new_name|
              keys = @models.build_model(model_name).build_template.build_attributes.attribute_keys
              if keys.include?(new_name)
                @errors << "Model #{model_name} trying to rename #{old_name} to #{new_name}, but a different #{new_name} already exists."
              end

              if !keys.include?(old_name)
                @errors << "Model #{model_name} trying to rename #{old_name} to #{new_name}, but #{old_name} does not exist."
              end
            end
          end
        end
      end

      class AddModelValidator < ValidatorBase
        attr_reader :model

        def initialize(models = Models.new, model_name = 'project')
          super()
          @model = models.build_model(model_name)
          @models = models
        end

        def parent_reciprocal_attribute
          @models.find_reciprocal(model: @model, link_attribute_name: @model.template.parent)
        end

        def name
          model&.name
        end

        def raw
          model&.template&.raw
        end

        def is_project?
          name == 'project'
        end

        def validate
          @errors << "Model name #{name} must be snake_case and can only consist of letters and \"_\"." unless name =~ name_regex_no_numbers

          if parent_reciprocal_attribute&.attribute_type != Etna::Clients::Magma::ParentLinkType::TABLE
            check_key("model #{name}", raw, 'identifier')
          end

          validate_links
          validate_attributes
        end

        def validate_attributes
          model.template.attributes.attribute_keys.each do |attribute_name|
            attribute = model.template.attributes.attribute(attribute_name)

            reciprocal = @models.find_reciprocal(model: model, attribute: attribute)
            if attribute_name == model.template.identifier && reciprocal&.attribute_type != AttributeType::TABLE
              attribute_types = [AttributeType::IDENTIFIER]
            elsif attribute_name == model.template.parent
              attribute_types = [AttributeType::PARENT]
            elsif reciprocal&.attribute_type == AttributeType::PARENT
              attribute_types = AttributeValidator.valid_parent_link_attribute_types
            else
              attribute_types = AttributeValidator.valid_add_row_attribute_types
            end

            attribute_validator = AttributeValidator.new(attribute, attribute_types, @models)
            attribute_validator.validate
            @errors += attribute_validator.errors
          end
        end

        def validate_links
          if !is_project? && !@models.model_keys.include?(@model.template.parent)
            @errors << "Parent model \"#{@model.template.parent}\" for #{name} does not exist in project."
          end

          if !is_project? && parent_reciprocal_attribute.nil?
            @errors << "Parent link attributes not defined for model #{name}."
          end

          link_attributes.each do |attribute|
            check_key("attribute #{attribute.attribute_name}", attribute.raw, 'link_model_name')

            # if attribute.attribute_name != attribute.link_model_name
            #   @errors << "Linked model, \"#{attribute.link_model_name}\", does not match attribute #{attribute.attribute_name}, link attribute names must match the model name."
            # end

            unless @models.model_keys.include?(attribute.link_model_name)
              @errors << "Linked model, \"#{attribute.link_model_name}\", on attribute #{attribute.attribute_name} does not exist!"
            end

            reciprocal = @models.find_reciprocal(model: @model, attribute: attribute)
            if reciprocal.nil?
              @errors << "Linked model, \"#{attribute.link_model_name}\", on attribute #{attribute.attribute_name} does not have a reciprocal link defined."
            end
          end
        end

        def link_attributes
          @model.template.attributes.all.select do |attribute|
            attribute.link_model_name
          end
        end
      end

      class AttributeValidator < ValidatorBase
        attr_reader :attribute

        def initialize(attribute, valid_attribute_types, project_models)
          super()
          @attribute = attribute
          @valid_attribute_types = valid_attribute_types
          @valid_validation_types = self.class.valid_validation_types
          @project_models = project_models
        end

        def self.valid_add_row_attribute_types
          Etna::Clients::Magma::AttributeType.entries.reject { |a|
            a == Etna::Clients::Magma::AttributeType::CHILD ||
                a == Etna::Clients::Magma::AttributeType::IDENTIFIER ||
                a == Etna::Clients::Magma::AttributeType::PARENT
          }.sort
        end

        def self.valid_parent_link_attribute_types
          [
              AttributeType::COLLECTION,
              AttributeType::TABLE,
              AttributeType::CHILD,
          ]
        end

        def self.valid_update_attribute_types
          Etna::Clients::Magma::AttributeType.entries.sort
        end

        def self.valid_validation_types
          Etna::Clients::Magma::AttributeValidationType.entries.sort
        end

        def validate
          validate_basic_attribute_data
          validate_attribute_validation
        end

        def validate_basic_attribute_data
          check_valid_name_with_numbers('Attribute', attribute.attribute_name)
          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'attribute_type')

          if attribute.link_model_name && ![AttributeType::TABLE, AttributeType::LINK, AttributeType::COLLECTION, AttributeType::PARENT, AttributeType::CHILD].include?(attribute.attribute_type)
            @errors << "attribute #{attribute.attribute_name} has link_model_name set, but has attribute_type #{attribute.attribute_type}"
          end

          if attribute.link_model_name && !@project_models.model_keys.include?(attribute.link_model_name)
            @errors << "attribute #{attribute.attribute_name} has link_model_name value of #{attribute.link_model_name}, but a model by that name does not exist."
          end

          check_in_set("attribute #{attribute.attribute_name}", attribute.raw, 'attribute_type', @valid_attribute_types)
        end

        def validate_attribute_validation
          return unless attribute.validation
          check_key("attribute #{attribute.attribute_name}, validation", attribute.validation, 'type')
          check_key("attribute #{attribute.attribute_name}, validation", attribute.validation, 'value')
          check_in_set("attribute #{attribute.attribute_name}, validation", attribute.validation, 'type', @valid_validation_types)
        end

        def is_link_attribute?
          attribute.attribute_type == Etna::Clients::Magma::AttributeType::LINK
        end

        def is_identifier_attribute?
          attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
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
          validator = AttributeValidator.new(@attribute, AttributeValidator.valid_add_row_attribute_types, project_models)
          validator.validate
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
