class Magma
  module Gnomon
    class Tokens
      include Enumerable

      def self.to_schema
        {
          gnomon_tokens: {
            type: "object",
            patternProperties: {
              "^[A-Z]+$": {
                "$ref": "#/definitions/gnomon_token"
              }
            }
          },
          gnomon_token: {
            type: "object",
            properties: {
              label: { type: "string" },
              values: {
                type: "object",
                patternProperties: {
                  ".*": { type: "string" }
                },
                minProperties: 1
              }
            },
            required: [ "values", "label" ]
          }
        }
      end

      attr_reader :errors

      def initialize(grammar)
        @grammar = grammar
        @errors = []
      end

      def valid?
        validations.each do |validation|
          send(validation)
        end

        @errors.empty?
      end

      def validations
        [
          :validate_single_separator
        ]
      end

      def config
        @grammar.config['tokens'] || {}
      end

      def each
        config.each do |token_name, token_values|
          yield token_name, token_values
        end
      end

      def match_value_by_description(token_name, match_token_name, match_token_value)
        return config[token_name]['values'].keys.find do |token_value|
          config[token_name]['values'][token_value] == config[match_token_name]['values'][match_token_value]
        end
      end

      def token_names
        config&.keys || []
      end

      def values(token_name)
        config[token_name]&.fetch('values')&.keys || []
      end

      def descriptions(token_name)
        config[token_name]&.fetch('values')&.values || []
      end

      def with_name(token_name)
        config[token_name]&.merge(
          name: token_name
        )
      end

      def valid_token_names(token_names)
        config.keys & token_names
      end

      def valid_tokens(token_names)
        valid_slice(token_names).values
      end

      def valid_values(token_names)
        valid_tokens(token_names).map do |token|
          token['values']
        end
      end

      def valid_descriptions(token_names)
        valid_tokens(token_names).map do |token|
          token['values'].values
        end
      end

      def valid_slice(token_names)
        config.slice( *valid_token_names(token_names) )
      end

      def is_separator?(token)
        Magma::Gnomon::Grammar::SEPARATOR_TOKENS.include?(token)
      end

      def is_numeric?(token)
        @grammar.numeric_increment == token || token =~ /_counter$/
      end

      private

      def validate_single_separator
        defined_separators = valid_token_names(Magma::Gnomon::Grammar::SEPARATOR_TOKENS)
        @errors << "No separator token defined!" unless defined_separators.size > 0
        @errors << "More than one separator token defined!" unless defined_separators.size < 2
        @errors << "More than one separator token defined!" unless defined_separators.all? {|sep| config[sep].fetch('values').size == 1}
      end
    end
  end
end
