module Redcap
  class Entity
    def self.to_schema
      {
        each: {
          type: "array",
          items: {
            oneOf: [
              { enum: ["record", "event", "repeat"] },
              { "$ref": "#/definitions/each_entity" },
            ]
          }
        },
        each_entity: {
          type: "object",
          properties: {
            record: { type: "string" },
            event: { type: "string" },
            repeat: { type: "string" },
            field: { type: "string" }
          },
          maxProperties: 1,
          minProperties: 1
        },
      }
    end

    class Record < Entity
      def key(eav)
        filtered_key(eav[:record])
      end
    end
    class Event < Entity
      def key(eav)
        filtered_key(eav[:redcap_event_name])
      end
    end
    class Repeat < Entity
      def key(eav)
        filtered_key([ eav[:redcap_repeat_instrument], eav[:redcap_repeat_instance] ], eav[:redcap_repeat_instrument])
      end
    end
    class Field < Entity
      def key(eav)
        if eav[:field_name] =~ @filter
          eav[:value]
        else
          :field
        end
      end

      def flat_key?
        false
      end
    end

    def self.create(entity)
      entity_name, filter = entity.is_a?(Hash) ? entity.first : [ entity, nil ]

      entity_class_name = entity_name.to_s.split('_').map(&:capitalize).join.to_sym

      raise "No such entity #{entity_name}" unless self.const_defined?(entity_class_name)

      self.const_get(entity_class_name).new(filter)
    end

    def initialize(filter)
      @filter = Regexp.new(filter) if filter
    end

    def filtered_key(key, value=nil)
      key unless @filter && (value || key) !~ @filter
    end

    def flat_key?
      true
    end
  end
end
