require_relative 'answer_collection_base'
require_relative 'answer_tuple'

class Magma
  class AnswerTupleArray < Magma::AnswerCollectionBase

    def self.from_raw_answer_tuples(raw_answer_tuples)
      Magma::AnswerTupleArray.new(
        raw_answer_tuples.map do |raw_answer_tuple|
          if raw_answer_tuple.is_a?(Magma::AnswerBase)
            raw_answer_tuple
          else
            inner_data = raw_answer_tuple.last

            if (!inner_data.is_a?(Magma::AnswerBase))
              if inner_data.is_a?(Array) # May break on file_collection queries ...
                inner_data = Magma::AnswerTupleArray.from_raw_answer_tuples(
                  inner_data
                )
              else
                inner_data = Magma::Answer.new(inner_data)
              end
            end

            Magma::AnswerTuple.new(
              raw_answer_tuple.first,
              inner_data
            )
          end
        end
      )
    end

    # def self.answer_tuple_array?(value)
    #   value.is_a?(Array) && value.all? do |val|
    #     Magma::AnswerTuple.answer_tuple?(val)
    #   end
    # end

    # def aggregated_values(data_is_collection)
    #   aggregate_nested_values(data_is_collection)
    # end

    # private

    # def aggregate_nested_values(data_is_collection)
    #   [].tap do |result|
    #     require 'pry'
    #     binding.pry
    #     queue = @raw_data.dup

    #     while !queue.empty?
    #       next_data = queue.shift
    #       if next_data.is_a?(Magma::AnswerTupleArray)
    #         queue += next_data.data
    #       elsif next_data.is_a?(Magma::AnswerTuple)
    #         if data_is_collection && !next_data.data.is_a?(Array)
    #           result = result.concat(next_data.data)
    #         else
    #           queue << next_data.data
    #         end
    #       else
    #         result << next_data.data
    #       end
    #     end
    #   end.flatten.compact
    # end
  end
end
