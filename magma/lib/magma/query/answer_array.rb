require_relative 'answer_collection_base'

class Magma
  class AnswerArray < Magma::AnswerCollectionBase
    def self.from_raw_answers(raw_answers)
      Magma::AnswerArray.new(
        raw_answers.map do |raw_answer|
          Magma::Answer.new(raw_answer)
        end
      )
    end
  end
end