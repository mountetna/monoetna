require_relative 'answer_base'

class Magma
  class AnswerTuple < Magma::AnswerBase

    # def self.answer_tuple_format?(value)
    #   # This is a weak / brittle area...
    #   (value.is_a?(Sequel::Postgres::JSONArray) || value.is_a?(Array)) &&
    #   value.length == 2 &&
    #   (value.first.is_a?(String) || value.first.is_a?(Integer)) # tables have integer indices
    # end

    attr_reader :identifier

    def initialize(identifier, data)
      @identifier = identifier
      @raw_data = data
      # parse_data
    end

    # def final_value
    #   @value
    # end

    def data
      @raw_data
    end

    def aggregated_values(data_is_collection)
      aggregate_nested_values(data_is_collection)
    end

    private

    def aggregate_nested_values(data_is_collection)
      [].tap do |result|
        queue = @raw_data.dup

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
    # private

    # def parse_data
    #   queue = [@raw_data]
    #   while queue.length > 0
    #     next_data = queue.shift
    #     if Magma::AnswerTuple.answer_tuple?(next_data)
    #       queue << next_data.last
    #     else
    #       @value = next_data
    #       break
    #     end
    #   end
    # end
  end
end