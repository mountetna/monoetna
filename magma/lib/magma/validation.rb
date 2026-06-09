require_relative 'validation/model'

class Magma
  class Validation
    attr_reader :grammar

    def initialize(grammar)
      @models = {}
      @grammar = grammar
    end

    def validate(model, record_name, record)
      model_validation(model).validate(record_name, record) do |error|
        yield error
      end
    end

    private

    def model_validation model
      @models[model] ||= Magma::Validation::Model.new(model, self)
    end
  end
end
