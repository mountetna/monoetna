module Redcap
  class Value
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
        return field_value_array(redcap_record).join(@config[:combine] || ' ')
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
      if @config.has_key?(:equals)
        return value == @config[:equals]
      end

      if @config.has_key?(:match)
        return value =~ @config[:match]
      end

      if @config.has_key?(:in)
        return @config[:in].include?(value)
      end

      if @config.has_key?(:exists)
        return value && !value.empty?
      end

      true
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
