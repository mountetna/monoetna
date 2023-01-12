require_relative 'answer_base'

class Magma
  class Answer < Magma::AnswerBase

    def initialize(data)
      @raw_data = data
    end

    def data
      @raw_data
    end
  end
end