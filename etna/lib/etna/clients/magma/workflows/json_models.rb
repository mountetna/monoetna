require 'json'

module Etna
  module Clients
    class Magma
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
          errors << "Missing required key for #{label}: \"#{key}\"." if !raw.key?(key)
          errors << "Invalid empty #{key} for #{label}: \"#{raw[key]}\"." if raw.key?(key) && nil_or_empty?(raw[key])
        end

        def check_type(label, raw, key, valid_types)
          errors << "Invalid #{key} for #{label}: \"#{raw[key]}\".\nShould be one of #{valid_types}." if raw.key?(key) && !valid_types.include?(raw[key].strip)
        end
      end

      class JsonProject < JsonBase
        def initialize(filepath:)
          super()
          @raw = JSON.parse(File.read(filepath))
          validate
        end

        def name
          @raw['project_name']&.strip
        end

        def name_full
          @raw['project_name_full']&.strip
        end

        def models
          @models ||= @raw['models'].map { |model_name, model_def|
            JsonModel.new(model_name, model_def) }
        end

        def fetch_model(name)
          models.select { |model| model.name == name }.first
        end

        def models_by_parent
          @models_by_parent ||= models.group_by { |model| model.parent_model_name }
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
            #   ... trying to make pagination consistent.
            sorted_models += child_models.sort { |m1, m2|
              m1.name <=> m2.name }
          end
          sorted_models
        end

        def get_magma_models
          # Convert models to a Magma Model + set of Magma Attributes in a
          #   template.
          magma_models = Etna::Clients::Magma::Models.new
          models.each { |model|
            model_builder = magma_models.build_model(model.name)
            model.to_magma_model(model_builder)

            # Because the input JSON format doesn't specify the child
            #   or reciprocal link attributes, we'll need to add those
            #   in manually.
            add_child_attributes(model_builder, model)
          }

          add_reciprocal_link_attributes(magma_models)
          magma_models
        end

        def add_child_attributes(builder, model)
          return unless models_by_parent[model.name]

          models_by_parent[model.name].each { |child_model|
            attribute_builder = builder.build_template.build_attributes
            add_child_attribute(attribute_builder, child_model)
          }
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
          all_link_attributes { |model, attribute|
            reciprocal_model = fetch_model(attribute.link_model_name)
            model_builder = magma_models.build_model(reciprocal_model.name)
            attribute_builder = model_builder.build_template.build_attributes
            add_reciprocal_link_attribute(attribute_builder, model)
          }
        end

        def add_reciprocal_link_attribute(builder, model)
          builder.build_attribute(model.name).tap do |attribute|
            attribute.attribute_name = model.name
            attribute.name = model.name
            attribute.attribute_type = Etna::Clients::Magma::AttributeType::COLLECTION
            attribute.display_name = model.prettified_name
            attribute.desc = model.prettified_name
            attribute.link_model_name = model.name
          end
        end

        def validate
          validate_project_names
          validate_models
          validate_model_links
        end

        def validate_project_names
          check_key('root project', @raw, 'project_name')
          check_key('root project', @raw, 'project_name_full')
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

          all_link_attributes { |model, attribute|
            check_key("model #{model.name}, attribute #{attribute.name}", attribute.raw, 'link_model_name')

            # Check that the linked model exists.
            errors << "Model \"#{model.name}\" already belongs to parent model \"#{model.parent_model_name}\". Remove attribute \"#{attribute.name}\"." if attribute.link_model_name == model.parent_model_name

            errors << "Linked model, \"#{attribute.link_model_name}\", on attribute #{attribute.name} of model #{model.name} does not exist!" if !model_names.include?(attribute.link_model_name)
          }
        end

        def all_link_attributes
          models.map { |model|
            model.attributes.map { |attribute|
              yield [model, attribute] if attribute.type == Etna::Clients::Magma::AttributeType::LINK
            }
          }
        end
      end

      class JsonModel < JsonBase
        attr_reader :name
        def initialize(model_name, raw)
          super()
          @name = model_name.strip
          @raw = raw

          @valid_parent_link_types = [
            Etna::Clients::Magma::ParentLinkType::CHILD,
            Etna::Clients::Magma::ParentLinkType::COLLECTION,
            Etna::Clients::Magma::ParentLinkType::TABLE
          ]
          validate
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
          @attributes ||= @raw.key?('attributes') ? @raw['attributes'].map { |attribute_name, attribute_def|
            JsonAttribute.new(self, attribute_name, attribute_def) } : []
        end

        def validate
          validate_add_model_data
          validate_attributes
        end

        def validate_add_model_data
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
          attributes.each { |attribute|
            @errors += attribute.errors unless attribute.valid?
          }
        end

        def to_magma_model(builder)
          template_builder = builder.build_template
          template_builder.identifier = identifier
          template_builder.parent = parent_model_name
          attribute_builder = template_builder.build_attributes
          attributes.each { |attribute|
            attribute.to_magma_model(attribute_builder)
          }

          # Because the input JSON format doesn't specify the parent
          #   attribute, we'll need to add those
          #   in manually.
          add_parent_attribute(attribute_builder) if !is_project?
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

      class JsonAttribute < JsonBase
        attr_reader :name, :model, :raw
        def initialize(model, attribute_name, raw)
          super()
          @model = model
          @name = attribute_name.strip
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

          validate
        end

        def type
          @raw['attribute_type']&.strip
        end

        def link_model_name
          @raw['link_model_name']&.strip
        end

        def validate
          validate_add_attribute_data
        end

        def validate_add_attribute_data
          check_key("model #{model.name}, attribute #{name}", @raw, 'attribute_name')

          if model.identifier != name
            check_key("model #{model.name}, attribute #{name}", @raw, 'attribute_type')

            # The following two could be calculated or left blank?
            # But it would be nice in the final UI to have them be more informative, so
            #   we'll enforce here.
            check_key("model #{model.name}, attribute #{name}", @raw, 'display_name')
            check_key("model #{model.name}, attribute #{name}", @raw, 'desc')

            check_type("model #{model.name}, attribute #{name}", @raw, 'attribute_type', @valid_attribute_types)
          end

          if @raw.key?('validation')
            check_key("model #{model.name}, attribute #{name}, validation", @raw['validation'], 'type')
            check_key("model #{model.name}, attribute #{name}, validation", @raw['validation'], 'value')
            check_type("model #{model.name}, attribute #{name}, validation", @raw['validation'], 'type', @valid_validation_types)
          end
        end

        def to_magma_model(builder)
          builder.build_attribute(name).tap do |attribute|
            attribute.attribute_name = name
            attribute.name = name

            if model.identifier != name
              attribute.display_name = @raw['display_name'].strip
              attribute.desc = @raw['desc'].strip
              attribute.attribute_type = type
            else
              attribute.attribute_type = Etna::Clients::Magma::AttributeType::IDENTIFIER
            end

            attribute.hidden = @raw['hidden'] if @raw.key?('hidden')
            attribute.read_only = @raw['read_only'] if @raw.key?('read_only')
            attribute.validation = @raw['validation'] if @raw.key?('validation')
            attribute.restricted = @raw['restricted'] if @raw.key?('restricted')
            attribute.link_model_name = link_model_name if link_model_name
          end
        end
      end
    end
  end
end
