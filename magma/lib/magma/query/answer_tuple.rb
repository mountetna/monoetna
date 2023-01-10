class Magma
  class AnswerTuple

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
        if answer_tuple?(next_data)
          queue << next_data.last
        else
          @value = next_data
          break
        end
      end
    end

    def answer_tuple?(value)
      (value.is_a?(Sequel::Postgres::JSONArray) || value.is_a?(Array)) &&
      value.length == 2 &&
      value.last.is_a?(Array)
    end
  end
end