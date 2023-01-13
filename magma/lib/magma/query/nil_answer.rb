require_relative 'simple_answer_base'

class Magma
  class NilAnswer < Magma::SimpleAnswerBase
    def to_s
      nil
    end
  end
end