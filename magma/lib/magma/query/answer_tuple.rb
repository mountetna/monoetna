class Magma
  class AnswerTuple

    def self.answer_tuple?(value)
      (value.is_a?(Sequel::Postgres::JSONArray) || value.is_a?(Array)) &&
      value.length == 2 &&
      value.first.is_a?(String)
    end

    def initialize(data)
      @raw_data = data
      parse_data
    end

    def final_value
      @value
    end

    private

    def parse_data
      queue = [@raw_data]
      while queue.length > 0
        next_data = queue.shift
        if Magma::AnswerTuple.answer_tuple?(next_data)
          queue << next_data.last
        else
          @value = next_data
          break
        end
      end
    end
  end
end