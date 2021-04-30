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
      else
        raise "Invalid value specified in value config: #{@config}"
      end
    end

    def valid?(record)
      if @config.has_key?(:equals)
        return field_value(record) == @config[:equals]
      end

      if @config.has_key?(:match)
        return field_value(record) =~ @config[:match]
      end

      if @config.has_key?(:in)
        return @config[:in].include?(field_value(record))
      end

      if @config.has_key?(:exists)
        return field_value(record) && !field_value(record).empty?
      end

      true
    end

    private

    def field_value(record)
      record[ field_name&.to_sym || @attribute_name.to_s ]
    end
  end
end
