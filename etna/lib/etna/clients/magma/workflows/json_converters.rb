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
