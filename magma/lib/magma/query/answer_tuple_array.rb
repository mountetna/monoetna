require_relative 'answer_tuple'

class Magma
  class AnswerTupleArray

    def self.answer_tuple_array?(value)
      value.is_a?(Array) && value.all? do |val|
        Magma::AnswerTuple.answer_tuple?(val)
      end
    end

    def initialize(raw_data)
      @raw_data = raw_data
    end

    def aggregated_values(data_is_collection)
      return [] unless Magma::AnswerTupleArray.answer_tuple_array?(@raw_data)

      aggregate_nested_values(data_is_collection)
    end

    private

    def aggregate_nested_values(data_is_collection)
      [].tap do |result|
        queue = @raw_data.dup

        while !queue.empty?
          next_data = queue.shift
          if Magma::AnswerTupleArray.answer_tuple_array?(next_data)
            queue += next_data
          elsif Magma::AnswerTuple.answer_tuple?(next_data)
            if data_is_collection && !next_data.last.is_a?(Array)
              result = result.concat(next_data)
            else
              queue << next_data.last
            end
          else
            result << next_data
          end
        end
      end.flatten.compact
    end
  end
end
