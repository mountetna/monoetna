require_relative 'answer_base'

class Magma
  class AnswerAggregationBase < Magma::AnswerBase

    def aggregated_values(data_is_collection=false)
      aggregate_nested_values(data_is_collection)
    end

    private

    def aggregate_nested_values(data_is_collection)
      [].tap do |result|
        queue = @raw_data.dup

        queue = [ queue ] unless queue.is_a?(Array)

        while !queue.empty?
          next_data = queue.shift
          if next_data.is_a?(Magma::AnswerTupleArray)
            queue += next_data.data
          elsif next_data.is_a?(Magma::AnswerTuple)
            if data_is_collection && !next_data.data.is_a?(Array)
              result = result.concat(next_data.data)
            elsif next_data.data.is_a?(Array)
              queue += next_data.data
            else
              queue << next_data.data
            end
          else
            result << next_data.data
          end
        end
      end.flatten.compact
    end
  end
end