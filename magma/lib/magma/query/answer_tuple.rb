require_relative 'answer_aggregation_base'

class Magma
  class AnswerTuple < Magma::AnswerAggregationBase
    attr_reader :identifier

    def initialize(identifier, data)
      @identifier = identifier
      @raw_data = data
    end

    def data
      @raw_data
    end

    def to_json(opts=nil)
      data_value = data.to_json(opts)
      begin
        data_value = JSON.parse(data_value)
      rescue TypeError, JSON::ParserError => e
      end

      [ @identifier, data_value ].to_json
    end
  end
end