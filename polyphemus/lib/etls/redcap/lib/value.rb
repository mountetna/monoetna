module Redcap
  class Value
    def self.to_schema
      {
        attribute_value: {
          type: "object",
          properties: {
            redcap_field: { type: "string" },
            value: { enum: [ "text", "value", "label", "note", "select_choice", "combine", "age", "none" ] },
            text: { type: "string" },
            combine: { type: "string" },
            equals: { type: "string" },
            match: { type: "string" },
            in: { type: "array", items: { type: "string" } },
            exists: { type: "boolean" }
          },
          additionalProperties: false,
          required: [ "value" ]
        }
      }
    end

    def initialize(att_name, config, template)
      @attribute_name = att_name
      @template = template
      @config = config.is_a?(String) ? { redcap_field: config, value: 'value' } : config
    end

    def field_name
      @config[:redcap_field]
    end

    def none?
      @config[:value] == "none"
    end

    def to_value(redcap_record, id)
      case @config[:value]
      when "text"
        return @config[:text]
      when "value"
        return field_value(redcap_record)
      when "label"
        return @template.field_label(field_name)
      when "note"
        return @template.field_note(field_name)
      when "select_choice"
        return @template.select_choice(field_name, field_value(redcap_record))
      when "combine"
        return field_value_array(redcap_record)&.join(@config[:combine] || ' ')
      when "age"
        raw_value = field_value(redcap_record)
        return [ 89, raw_value.to_i ].min.to_s if raw_value
      else
        raise "Invalid value specified in value config: #{@config}"
      end
    end

    def valid_redcap?(redcap_record)
      valid?(field_value(redcap_record))
    end

    def valid_magma?(magma_record)
      valid?(magma_value(magma_record))
    end

    def valid?(value)
      # At some point we would ideally deprecate this filtering ability in the value,
      #   but that would require migrating current scripts. This at least consolidates
      #   the validity logic.
      @filter = Redcap::Filter.new(@config)

      @filter.valid?(value)
    end

    private

    def field_value_array(redcap_record)
      redcap_record[ field_name&.to_sym ]
    end

    def magma_value(magma_record)
      magma_record[ @attribute_name.to_s ]
    end

    def field_value(redcap_record)
      field_value_array(redcap_record)&.first
    end
  end
end
