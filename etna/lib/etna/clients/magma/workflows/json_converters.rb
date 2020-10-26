require 'json'

module Etna
  module Clients
    class Magma
      class ConverterBase
        def self.convert_attribute_user_json_to_magma_json(user_model_json)
          json_attribute = Marshal.load(Marshal.dump(user_model_json['attributes']))
          return unless json_attribute

          json_attribute.keys.each do |attribute_name|
            json_attribute[attribute_name]['attribute_type'] = Etna::Clients::Magma::AttributeType::IDENTIFIER if user_model_json['identifier'] == attribute_name
          end
          json_attribute
        end

        def self.convert_model_user_json_to_magma_json(model_name, user_json)
          json_model = Marshal.load(Marshal.dump(user_json))
          json_model['template'] = {
            'name' => model_name,
            'identifier' => user_json['identifier'],
            'parent' => user_json['parent_model_name'],
            'attributes' => convert_attribute_user_json_to_magma_json(json_model)
          }
          json_model
        end

        def self.convert_project_user_json_to_magma_json(user_json)
          magma_models_json = {}
          user_json['models'].keys.each do |model_name|
            magma_models_json[model_name] = convert_model_user_json_to_magma_json(
              model_name,
              user_json['models'][model_name])
          end
          user_json['models'] = magma_models_json
          user_json
        end

        def prettify(name)
          name.split('_').map(&:capitalize).join(' ')
        end
      end

      class ProjectConverter < ConverterBase
        # Add in missing attributes from user JSON -> Magma JSON
        attr_reader :project
        def initialize(project)
          @project = project
        end

        def models_by_parent
          @models_by_parent ||= project.models.all.group_by { |model| model.raw['parent_model_name'] }
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

        def convert!
          project.models.all.each do |model|
            # Because the input JSON format doesn't specify the child
            #   or reciprocal link attributes, we'll need to add those
            #   in manually.
            add_child_attributes(model)

            model_converter = ModelConverter.new(model)
            model_converter.convert!
          end

          add_reciprocal_link_attributes
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
            attribute.attribute_type = child_model.raw['parent_link_type']
            attribute.display_name = prettify(child_model.name)
            attribute.desc = prettify(child_model.name)
          end
        end

        def add_reciprocal_link_attributes
          all_link_attributes do |model, attribute|
            reciprocal_model = project.models.model(attribute.link_model_name)
            model_converter = ModelConverter.new(reciprocal_model)
            model_converter.add_reciprocal_link_attribute(model)
          end
        end

        def all_link_attributes
          project.models.all.map do |model|
            model.template.attributes.all.select do |attribute|
              attribute.attribute_type == Etna::Clients::Magma::AttributeType::LINK
            end.each do |attribute|
              yield [model, attribute]
            end
          end
        end
      end

      class ModelConverter < ConverterBase
        attr_reader :model
        def initialize(model)
          @model = model
        end

        def is_project?
          model.name == 'project'
        end

        def parent_model_name
          model.raw['parent_model_name']
        end

        def identifier
          model.raw['identifier']
        end

        def add_reciprocal_link_attribute(reciprocal_model)
          attribute_builder = model.build_template.build_attributes
          attribute_builder.build_attribute(reciprocal_model.name).tap do |attribute|
            attribute.attribute_name = reciprocal_model.name
            attribute.name = reciprocal_model.name
            attribute.attribute_type = Etna::Clients::Magma::AttributeType::COLLECTION
            attribute.display_name = prettify(reciprocal_model.name)
            attribute.desc = prettify(reciprocal_model.name)
            attribute.link_model_name = reciprocal_model.name
          end
        end

        def convert!
          template_builder = model.build_template
          template_builder.identifier = identifier
          template_builder.parent = parent_model_name
          attribute_builder = template_builder.build_attributes

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
            attribute.display_name = prettify(parent_model_name)
            attribute.desc = prettify(parent_model_name)
          end
        end
      end

      class AttributeActionsConverter < ConverterBase
        attr_reader :actions
        def initialize(actions)
          @actions = JSON.parse(actions.to_json, symbolize_names: true)
        end

        def camelize(action_name)
          action_name.split('_').map(&:capitalize).join('')
        end

        def clazz_name(action_name)
          "Etna::Clients::Magma::#{camelize(action_name)}Action"
        end

        def convert
          actions.map do |action_json|
            # We use desc and attribute_type to be consistent
            #   with the other JSON actions...but the
            #   Magma Model takes type and description.
            if action_json[:action_name] == 'add_attribute'
              action_json[:type] = action_json.delete(:attribute_type)
              action_json[:description] = action_json.delete(:desc)
            elsif action_json[:action_name] == 'add_link'
              action_json[:links].first[:type] = action_json[:links].first.delete(:attribute_type)
              action_json[:links].last[:type] = action_json[:links].last.delete(:attribute_type)
            end

            clazz = Object.const_get(clazz_name(action_json[:action_name]))
            clazz.new(**action_json)
          rescue ArgumentError => e
            modified_message = "Exception while parsing #{action_json}.\n" + e.message
            raise ArgumentError.new(modified_message)
          end
        end
      end
    end
  end
end
