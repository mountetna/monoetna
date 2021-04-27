module Redcap
  class Value
    def initialize(config, template)
      @template = template
      @config = config.is_a?(String) ? { redcap_field: config, value: 'value' } : config
    end

    def field_name
      @config[:redcap_field]
    end

    def to_value(redcap_record, id)
      case @config[:value]
      when "text"
        return @config[:text]
      when "value"
        return field_value(redcap_record)
      when "label"
        return @template.field_label(@config[:redcap_field])
      when "note"
        return @template.field_note(@config[:redcap_field])
      else
        raise "Invalid value specified in value config: #{@config}"
      end
    end

    def valid?(redcap_record)
      return field_value(redcap_record) == @config[:equals] if @config.has_key?(:equals)
      return @config[:in].include?(field_value(redcap_record)) if @config.has_key?(:in)
      return field_value(redcap_record) && !field_value(redcap_record).empty? if @config.has_key?(:exists)

      true
    end

    private

    def field_value(redcap_record)
      redcap_record[ @config[:redcap_field].to_sym ]
    end
  end
end
