require_relative 'answer_base'

class Magma
  class AnswerCollectionBase < Magma::AnswerBase
    def initialize(raw_data)
      @raw_data = raw_data
    end

    def data
      @raw_data
    end

    def [](key)
      @raw_data[key]
    end

    def map(&block)
      @raw_data.map(&block)
    end

  end
end