require_relative 'answer_aggregation_base'

class Magma
  class AnswerCollectionBase < Magma::AnswerAggregationBase
    def self.array_of_answers?(data)
      data.is_a?(Array) && data.all? do |datum|
        datum.is_a?(Magma::AnswerBase)
      end
    end

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

    def to_json(opts=nil)
      @raw_data.map do |d|
        data_value = d.to_json(opts)

        begin
          data_value = JSON.parse(data_value)
        rescue TypeError, JSON::ParserError => e
        end

        data_value
      end.to_json
    end

    def to_h
      JSON.parse(to_json).to_h
    end
  end
end