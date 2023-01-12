class Magma
  module WithAggregatedDataModule
    def nested_reduce_and_apply(nested_array_data, level, &transform)
      return [Magma::Answer.new(transform.call(nested_array_data))] unless nested_array_data.respond_to?(:reduce)

      nested_array_data&.reduce([]) do |result, value|
        if include_identifier?(value, level)
          # is an tuple, [identifier, data]
          identifier = value.first
          raw_data = value.last

          if (raw_data.nil?)
            data = [Magma::NilAnswer.new]
          elsif !(raw_data.is_a?(Array) || raw_data.is_a?(Sequel::Postgres::JSONArray))
            data = [Magma::Answer.new(transform.call(raw_data))]
          else
            data = nested_reduce_and_apply(raw_data, level + 1, &transform)
          end

          result << Magma::AnswerTuple.new(identifier, data)
        elsif skip_identifier?(value, level)
          # Skip the intervening identifiers
          raw_data = value.last
          if (raw_data.nil?)
            result << Magma::NilAnswer.new
          else
            result = result.concat(nested_reduce_and_apply(raw_data, level + 1, &transform))
          end
        else
          require 'pry'
          binding.pry
          result << Magma::Answer.new(transform.call(value))
        end
      end
    end

    def include_identifier?(value, level)
      value.is_a?(Array) && 0 == level
    end

    def skip_identifier?(value, level)
      value.is_a?(Array) && level > 0
    end
  end
end