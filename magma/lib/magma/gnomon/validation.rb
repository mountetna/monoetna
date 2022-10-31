class Magma
  module Gnomon
    class Validation
      attr_reader :errors

      def initialize(grammar)
        @grammar = grammar
        @errors = []
      end

      def valid?
        validations.each do |validation|
          send(validation)
        end

        errors.empty?
      end

      private

      def validations
        [
          :validate_schema,
          :validate_rules,
          :validate_tokens,
          :validate_synonyms
        ]
      end

      def validate_schema
        schema = JSONSchemer.schema(
          JSON.parse(Magma::Gnomon::Grammar.to_schema.to_json)
        )

        schema_errors = schema.validate(JSON.parse(@grammar.config.to_json))

        @errors += schema_errors.map do |error|
          JSONSchemer::Errors.pretty(error)
        end
      end

      def validate_rules
        @errors += @grammar.rules.errors unless @grammar.rules.valid?
      end

      def validate_tokens
        @errors += @grammar.tokens.errors unless @grammar.tokens.valid?
      end

      def validate_synonyms
        @errors += @grammar.synonyms.errors unless @grammar.synonyms.valid?
      end
    end
  end
end
