require_relative 'answer_collection_base'
require_relative 'answer_tuple'

class Magma
  class AnswerTupleArray < Magma::AnswerCollectionBase
    def self.from_raw_answer_tuples(raw_answer_tuples)
      Magma::AnswerTupleArray.new(
        raw_answer_tuples.map do |raw_answer_tuple|
          if (raw_answer_tuple.is_a?(Magma::AnswerBase) ||
              Magma::AnswerCollectionBase.array_of_answers?(raw_answer_tuple))
            raw_answer_tuple
          else
            identifier = raw_answer_tuple.first
            inner_data = raw_answer_tuple.last

            if (!inner_data.is_a?(Magma::AnswerBase))
              if inner_data.is_a?(Array)
                inner_data = Magma::AnswerTupleArray.from_raw_answer_tuples(
                  inner_data
                )
              else
                inner_data = Magma::Answer.new(inner_data)
              end
            end

            Magma::AnswerTuple.new(
              identifier,
              inner_data
            )
          end
        end
      )
    end
  end
end
