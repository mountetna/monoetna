module Redcap
  class Filter
    def self.to_schema
      {
        filter_value: {
          type: "object",
          properties: {
            redcap_field: { type: "string" },
            equals: { type: "string" },
            match: { type: "string" },
            in: { type: "array", items: { type: "string" } },
            exists: { type: "boolean" },
          },
          additionalProperties: false,
          required: [ "redcap_field" ],
          maxProperties: 2,
          minProperties: 2
        },
      }
    end

    def self.create(filter)
      Redcap::Filter.new(filter)
    end

    def initialize(config)
      @config = config
    end

    def field_name
      @config[:redcap_field]
    end

    def allow_redcap?(redcap_record)
      valid?(field_value(redcap_record))
    end

    def valid?(value)
      if @config.has_key?(:equals)
        return value == @config[:equals]
      end

      if @config.has_key?(:match)
        return value =~ Regexp.new(@config[:match])
      end

      if @config.has_key?(:in)
        return @config[:in].include?(value)
      end

      if @config.has_key?(:exists)
        return @config[:exists] ? (value && !value.empty?) : (value.nil? || value.empty?)
      end

      true
    end

    private

    def field_value_array(redcap_record)
      redcap_record[field_name&.to_sym]
    end

    def field_value(redcap_record)
      field_value_array(redcap_record)&.first
    end
  end
end
