require_relative 'answer_base'

class Magma
  class SimpleAnswerBase < Magma::AnswerBase
    def initialize(data=nil)
      @raw_data = data
    end

    def data
      @raw_data
    end

    def to_json(opts=nil)
      begin
        data.to_json
      rescue TypeError, JSON::ParserError => e
        data
      end
    end

    def to_s
      Magma::AnswerCollectionBase.array_of_answers?(data) ?
        data.map do |datum|
          datum.to_s
        end :
        data.to_s
    end
  end
end