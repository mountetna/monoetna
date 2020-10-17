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
          @errors << "Missing required key for #{label}: \"#{key}\"." if !raw.key?(key)
          @errors << "Invalid empty #{key} for #{label}: \"#{raw[key]}\"." if raw.key?(key) && nil_or_empty?(raw[key])
        end

        def check_type(label, raw, key, valid_types)
          @errors << "Invalid #{key} for #{label}: \"#{raw[key]}\".\nShould be one of #{valid_types}." if raw.key?(key) && !valid_types.include?(raw[key].strip)
        end

        def name_regex_with_numbers
          /\A[a-z][a-z0-9]*(_[a-z0-9]+)*\Z/
        end

        def name_regex_no_numbers
          /\A[a-z]*(_[a-z]+)*\Z/
        end

        def model_exists_in_project?(project_magma_models, model_name)
          !!project_magma_models.model(model_name)
        end

        def convert_model_user_json_to_magma_json(model_name, user_json)
          json_model = user_json.dup
          json_model['template'] = {
            'name' => model_name,
            'attributes' => json_model['attributes'].dup
          }
          json_model.delete('attributes')
          Etna::Clients::Magma::Model.new(json_model)
        end

        def convert_project_user_json_to_magma_json(user_json)
          magma_models_json = {}
          user_json['models'].keys.each do |model_name|
            magma_models_json[model_name] = convert_model_user_json_to_magma_json(
              model_name,
              user_json['models'][model_name])
          end
          user_json['models'] = magma_models_json
          user_json
        end
      end

      class ProjectValidator < ValidatorBase
        attr_reader :project
        def initialize(project)
          super
          @project = project
          validate
        end

        def models_by_parent
          @models_by_parent ||= models.all.group_by { |model| model.raw['parent_model_name'] }
        end

        def model_tree
          # Shallow -> deep
          sorted_models = []
          parent_model_names = [nil]

          loop do
            child_models = models_by_parent.values_at(*parent_model_names).flatten.compact

            break if child_models.length == 0

            parent_model_names = child_models.map { |model| model.name }

            # Sort by name within each depth level
            sorted_models += child_models.sort do |m1, m2|
              m1.name <=> m2.name
            end
          end
          sorted_models
        end

        def update_magma_models
          models.each do |model|
            # Because the input JSON format doesn't specify the child
            #   or reciprocal link attributes, we'll need to add those
            #   in manually.
            add_child_attributes(model_builder, model)
          end

          add_reciprocal_link_attributes(magma_models)
        end

        def add_child_attributes(model)
          return unless models_by_parent[model.name]

          models_by_parent[model.name].each do |child_model|
            attribute_builder = model.build_template.build_attributes
            add_child_attribute(attribute_builder, child_model)
          end
        end

        def add_child_attribute(builder, child_model)
          builder.build_attribute(child_model.name).tap do |attribute|
            attribute.attribute_name = child_model.name
            attribute.name = child_model.name
            attribute.attribute_type = child_model.parent_link_type
            attribute.display_name = child_model.prettified_name
            attribute.desc = child_model.prettified_name
          end
        end

        def add_reciprocal_link_attributes(magma_models)
          all_link_attributes do |model, attribute|
            reciprocal_model = project.models.model(attribute.link_model_name)
            model_builder = magma_models.build_model(reciprocal_model.name)
            reciprocal_model.add_reciprocal_link_attribute(model_builder, model)
          end
        end

        def validate
          validate_project_names
          validate_model_links
          validate_models
        end

        def validate_project_names
          check_key('root project', project.raw, 'project_name')
          check_key('root project', project.raw, 'project_name_full')
          name = project.raw['name']
          @errors << "Project name #{name} must be snake_case and cannot start with a number or \"pg_\"." unless name =~ name_regex_with_numbers && !name.start_with?('pg_')
        end

        def validate_models
          project.models.all.each do |model|
            validator = ModelValidator.new
            validator.validate_model(model)
            @errors += validator.errors unless validator.valid?
          end
        end

        def validate_model_links
          # Check all the attributes of type 'link', and make sure
          #   that the link_model_name exists in the project definition.
          model_names = project.models.all.map { |m| m.name }
          project.models.all.map do |model|
            validator = ModelValidator.new
            validator.validate_link_models(model_names)
            @errors += validator.errors unless validator.valid?
          end
        end

        def all_link_attributes
          models.map do |model|
            model.link_attributes do |attribute|
              yield [model, attribute]
            end
          end
        end
      end

      class ModelValidator < ValidatorBase
        attr_reader :name
        def initialize(model_name, raw)
          super()
          @name = model_name.strip
          @raw = raw

          @valid_parent_link_types = Etna::Clients::Magma::ParentLinkType.entries.sort # sort for prettier presentation

          validate
        end

        def self.linking_stub_from_name(model_name)
          # Used for linking to a link_model_name only! Fills out required fields with
          #   dummy information to pass validation.
          Etna::Clients::Magma::JsonModel.new(
            model_name,
            {
              "parent_model_name" => "stub",
              "parent_link_type" => "child",
              "identifier" => "none",
              "attributes" => {}
            }
          )
        end

        def parent_model_name
          @raw['parent_model_name']&.strip
        end

        def prettified_parent_model_name
          parent_model_name.split('_').map(&:capitalize).join(' ')
        end

        def parent_link_type
          @raw['parent_link_type']&.strip
        end

        def identifier
          @raw['identifier']&.strip
        end

        def prettified_name
          name.split('_').map(&:capitalize).join(' ')
        end

        def is_project?
          name == 'project'
        end

        def attributes
          @attributes ||= @raw.key?('attributes') ? @raw['attributes'].map do |attribute_name, attribute_def|
            attribute = Etna::Clients::Magma::Attribute.new(attribute_def)
            attribute.attribute_type = Etna::Clients::Magma::AttributeType::IDENTIFIER if identifier == attribute.attribute_name
            attribute
          end: []
        end

        def link_attributes
          attributes.map do |attribute|
            yield attribute if attribute.is_link_attribute?
          end
        end

        def add_reciprocal_link_attribute(builder, model)
          attribute_builder = builder.build_template.build_attributes
          attribute_builder.build_attribute(model.name).tap do |attribute|
            attribute.attribute_name = model.name
            attribute.name = model.name
            attribute.attribute_type = Etna::Clients::Magma::AttributeType::COLLECTION
            attribute.display_name = model.prettified_name
            attribute.desc = model.prettified_name
            attribute.link_model_name = model.name
          end
        end

        def validate
          validate_add_model_data
          validate_attributes
        end

        def validate_add_model_data
          @errors << "Model name #{name} must be snake_case and can only consist of letters and \"_\"." unless name =~ name_regex_no_numbers

          if !is_project?
            check_key("model #{name}", @raw, 'parent_model_name')
            check_key("model #{name}", @raw, 'parent_link_type')
            check_type("model #{name}", @raw, 'parent_link_type', @valid_parent_link_types)
          end
          if parent_link_type != Etna::Clients::Magma::ParentLinkType::TABLE
            check_key("model #{name}", @raw, 'identifier')
          end
        end

        def validate_attributes
          attributes.each do |attribute|
            attribute_validator = AttributeValidator.new
            attribute_validator.validate_add_attribute_data(attribute)
            @errors += attribute_validator.errors unless attribute_validator.valid?
          end
        end

        def validate_link_models(model_names)
          @errors << "Parent model #{parent_model_name} for #{name} does not exist in project.\nCurrent models are #{model_names}." if name != 'project' && !model_names.include?(parent_model_name)

          link_attributes do |attribute|
            attribute.validate_link_models(model_names)
            @errors = (@errors + attribute.errors).uniq # need to re-add attribute.errors because validate_link_models called at project level.

            @errors << "Model \"#{name}\" already belongs to parent model \"#{parent_model_name}\". Remove attribute \"#{attribute.name}\"." if attribute.link_model_name == parent_model_name
          end
        end

        def validate_existing_models(project_magma_models)
          @errors << "Model #{name} already exists in project!" if model_exists_in_project?(project_magma_models, name)
        end

        def to_magma_model(builder)
          template_builder = builder.build_template
          template_builder.identifier = identifier
          template_builder.parent = parent_model_name
          attribute_builder = template_builder.build_attributes
          attributes.each do |attribute|
            attribute.to_magma_model(attribute_builder)
          end

          # Because the input JSON format doesn't specify the parent
          #   attribute, we'll need to add those
          #   in manually.
          add_parent_attribute(attribute_builder) unless is_project?
        end

        def add_parent_attribute(builder)
          builder.build_attribute(parent_model_name).tap do |attribute|
            attribute.attribute_name = parent_model_name
            attribute.name = parent_model_name
            attribute.link_model_name = parent_model_name
            attribute.attribute_type = Etna::Clients::Magma::AttributeType::PARENT
            attribute.display_name = prettified_parent_model_name
            attribute.desc = prettified_parent_model_name
          end
        end
      end

      class AttributeValidator < ValidatorBase
        def initialize
          super
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

        def check_valid_name_with_numbers(label, proposed_name)
          @errors << "#{label} name #{proposed_name} must be snake_case and can only consist of letters, numbers, and \"_\"." unless proposed_name =~ name_regex_with_numbers
        end

        def validate_basic_attribute_data(attribute)
          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'attribute_type')

          # The following two could be calculated or left blank?
          # But it would be nice in the final UI to have them be more informative, so
          #   we'll enforce here.
          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'display_name')
          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'desc')

          check_type("attribute #{attribute.attribute_name}", attribute.raw, 'attribute_type', @valid_attribute_types)
        end

        def validate_attribute_validation(attribute)
          return unless raw.key?('validation')
          check_key("attribute #{attribute.attribute_name}, validation", attribute.validation, 'type')
          check_key("attribute #{attribute.attribute_name}, validation", attribute.validation, 'value')
          check_type("attribute #{attribute.attribute_name}, validation", attribute.validation, 'type', @valid_validation_types)
        end

        def is_link_attribute?(attribute)
          attribute.attribute_type == Etna::Clients::Magma::AttributeType::LINK
        end

        def is_identifier_attribute?(attribute)
          attribute.attribute_type == Etna::Clients::Magma::AttributeType::IDENTIFIER
        end

        def validate_link_models(attribute, model_names)
          return unless is_link_attribute?(attribute)

          check_key("attribute #{attribute.attribute_name}", attribute.raw, 'link_model_name')

          # Check that the linked model exists.
          @errors << "Linked model, \"#{attribute.link_model_name}\", on attribute #{attribute.attribute_name} does not exist!\nCurrent models are #{model_names}." unless model_names.include?(attribute.link_model_name)
        end

        def validate_add_attribute_data(attribute)
          check_valid_name_with_numbers('Attribute', attribute.attribute_name)

          validate_basic_attribute_data(attribute) unless is_identifier_attribute?(attribute)
          validate_attribute_validation(attribute)

          @errors << "Attribute name #{attribute.attribute_name} should match the link_model_name, \"#{attribute.link_model_name}\"." unless attribute.link_model_name && attribute.link_model_name == attribute.attribute_name
        end
      end

      # class JsonAttribute < ValidatorBase
      #   attr_reader :name, :model, :raw
      #   def initialize(model, attribute_name, raw)
      #     super()
      #     @model = model
      #     @name = attribute_name.strip
      #     @raw = raw

      #     validate
      #   end

      #   def type
      #     @raw['attribute_type']&.strip
      #   end

      #   def is_link_attribute?
      #     type == Etna::Clients::Magma::AttributeType::LINK
      #   end

      #   def link_model_name
      #     @raw['link_model_name']&.strip
      #   end

      #   def attribute_name
      #     name
      #   end

      #   def display_name
      #     @raw['display_name']&.strip
      #   end

      #   def desc
      #     @raw['desc']&.strip
      #   end

      #   def validate
      #     validate_add_attribute_data
      #   end

      #   def attribute_type
      #     return Etna::Clients::Magma::AttributeType::IDENTIFIER if model.identifier == name
      #     type
      #   end

      #   def hidden
      #     @raw['hidden']
      #   end

      #   def read_only
      #     @raw['read_only']
      #   end

      #   def validation
      #     @raw['validation']
      #   end

      #   def restricted
      #     @raw['restricted']
      #   end

      #   def format_hint
      #     @raw['format_hint']
      #   end

      #   def unique
      #     @raw['unique']
      #   end

      #   def validate_link_models(model_names)
      #     return unless is_link_attribute?

      #     check_key("model #{model.name}, attribute #{name}", raw, 'link_model_name')

      #     # Check that the linked model exists.
      #     @errors << "Linked model, \"#{link_model_name}\", on attribute #{name} of model #{model.name} does not exist!\nCurrent models are #{model_names}." unless model_names.include?(link_model_name)
      #   end

      #   def validate_add_attribute_data
      #     check_valid_name_with_numbers('Attribute', name)

      #     validate_basic_attribute_data(model.name, name, @raw) if model.identifier != name
      #     validate_attribute_validation(model.name, name, @raw)

      #     if link_model_name
      #       @errors << "Attribute name #{name} in model #{model.name} should match the link_model_name, \"#{link_model_name}\"." unless link_model_name == name
      #     end
      #   end

      #   def to_magma_model(builder)
      #     builder.build_attribute(name).tap do |attribute|
      #       Etna::Clients::Magma::Attribute.copy(self, attribute)
      #     end
      #   end
      # end

      class JsonAttributeAction < ValidatorBase
        def initialize(raw, project_models)
          @raw = raw
          @project_models = project_models
          validate
        end

        def name
          @raw['attribute_name']&.strip
        end

        def model_name
          @raw['model_name']&.strip
        end

        def validate
          raise "Subclasses must implement this method."
        end

        def exists_in_magma_model?(magma_model_name, attribute_name)
          !!project_models.model(magma_model_name).template.attributes.attribute(attribute_name)
        end

        def validate_model_exists(magma_model_name)
          @errors << "Model #{magma_model_name} does not exist in project." unless model_exists_in_project?(project_models, magma_model_name)
        end

        def check_exists_in_model(magma_model_name, attribute_name)
          @errors << "Attribute #{attribute_name} already exists in model #{magma_model_name}." if exists_in_magma_model?(magma_model_name, attribute_name)
        end
      end

      class JsonAddAttributeAction < JsonAttributeAction
        def initialize
          super
          @attribute = Etna::Clients::Magma::Attribute.new(@raw)
        end

        def validate
          validate_model_exists(model_name)
          validate_attribute_data
        end

        def validate_attribute_data
          validator = AttributeValidator.new
          validator.validate_basic_attribute_data(@attribute)
          validator.validate_attribute_validation(@attribute)
          @errors += validator.errors unless validator.valid?
        end
      end

      class JsonAddLinkAction < JsonAttributeAction
        def links
          @raw['links']
        end

        def source
          links.first
        end

        def dest
          links.last
        end

        def validate
          validate_links
          validate_both_models_exist
          validate_link_data
        end

        def validate_links
          check_key("action #{@raw}", @raw, "links")
          @errors << "Must include two link entries, each with \"model_name\", \"attribute_name\", and \"type\"." unless links.length == 2
          check_key("link #{source}", source, "model_name")
          check_key("link #{source}", source, "attribute_name")
          check_key("link #{source}", source, "type")
          check_key("link #{dest}", dest, "model_name")
          check_key("link #{dest}", dest, "attribute_name")
          check_key("link #{dest}", dest, "type")
        end

        def validate_both_models_exist
          validate_model_exists(source['model_name'])
          validate_model_exists(dest['model_name'])
        end

        def validate_link_data
          # Make sure the attribute names don't already exist in the models,
          #   and that the types are valid.
          link_types = Set.new([source['type'], dest['type']])
          expected_link_types = Set.new([
            Etna::Clients::Magma::AttributeTypes::LINK,
            Etna::Clients::Magma::AttributeTypes::COLLECTION
          ])
          @errors << "You must have one \"link\" and one \"collection\" type in the links." unless link_types == expected_link_types

          check_exists_in_model(source['model_name'], source['attribute_name'])
          check_exists_in_model(dest['model_name'], dest['attribute_name'])
        end
      end

      class JsonRenameAttributeAction < JsonAttributeAction
        def new_attribute_name
          @raw['new_attribute_name']&.strip
        end

        def validate
          validate_action
          validate_model_exists(model_name)
          validate_proposed_name
        end

        def validate_action
          check_key("action #{@raw}", @raw, "model_name")
          check_key("action #{@raw}", @raw, "attribute_name")
          check_key("action #{@raw}", @raw, "new_attribute_name")
        end

        def validate_proposed_name
          check_valid_name_with_numbers('New attribute', new_attribute_name)
          check_exists_in_model(model_name, new_attribute_name)
        end
      end

      class JsonUpdateAttributeAction < JsonAttributeAction
        def initialize
          super
          @attribute = Etna::Clients::Magma::Attribute.new(@raw)
        end

        def validate
          validate_model_exists(model_name)
          validate_attribute_data
        end

        def validate_attribute_data
          validator = AttributeValidator.new
          validator.validate_attribute_validation(@attribute)
          @errors += validator.errors unless validator.valid?
        end
      end

      class JsonAttributeActionFactory
        def camelize(action_name)
          action_name.split('_').map(&:capitalize).join('')
        end

        def clazz_name(action_name)
          "Json#{camelize(action_name)}Action"
        end

        def self.from_json(raw, project_models)
          clazz = Object.const_get(clazz_name(raw['action_name']))
          clazz.new(raw, project_models)
        end
      end
    end
  end
end
