require_relative 'simple_answer_base'

class Magma
  class Answer < Magma::SimpleAnswerBase

    def [](key)
      data[key]
    end
  end
end