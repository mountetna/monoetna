require_relative 'answer_base'

class Magma
  class NilAnswer < Magma::AnswerBase

    def data
      nil
    end
  end
end