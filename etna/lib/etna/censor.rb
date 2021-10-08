module Etna
  class Censor
    def initialize(log_redact_keys)
      @log_redact_keys = log_redact_keys
    end

    def redact_keys
      @log_redact_keys
    end

    def redact(key, value)
      # Redact any values for the supplied key values, so they
      #   don't appear in the logs.
      return compact(value) unless redact_keys

      if redact_keys.include?(key)
        return "*"
      elsif value.is_a?(Hash)
        redacted_value = value.map do |value_key, value_value|
          [value_key, redact(value_key, value_value)]
        end.to_h
        return redacted_value
      end

      return compact(value)
    end

    private

    def compact(value)
      value = value.to_s
      value = value[0..500] + "..." + value[-100..-1] if value.length > 600
      value
    end
  end
end
