class Magma
  module Gnomon
    class Grammar < Sequel::Model
      SYNONYM_SCHEMA= {
        gnomon_synonyms: {
          type: "array",
          items: {
            "$ref": "#/definitions/gnomon_synonym"
          }
        },
        gnomon_synonym: {
          type: "array",
          items: {
            type: "string",
            pattern: "^[A-Z]+$"
          }
        }
      }

      RULE_SCHEMA= {
        gnomon_rules: {
          type: "object",
          patternProperties: {
            "^\\w+$": {
              type: "string"
            }
          }
        }
      }

      TOKEN_SCHEMA= {
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

      def self.to_schema
        {
          "$schema": "http://json-schema.org/draft-07/schema#",
          title: "Gnomon Grammar",
          description: "This grammar defines a set of tokens and rules which can be used to match or create valid identifiers for a Mount Etna project.",
          definitions: [
            RULE_SCHEMA,
            TOKEN_SCHEMA,
            SYNONYM_SCHEMA,
          ].reduce(&:merge),
          type: "object",
          properties: {
            "tokens": {
              "$ref": "#/definitions/gnomon_tokens"
            },
            "synonyms": {
              "$ref": "#/definitions/gnomon_synonyms"
            },
            "rules": {
              "$ref": "#/definitions/gnomon_rules"
            }
          },
          required: [ "tokens", "rules" ],
          additionalProperties: false
        }
      end

      one_to_many :identifiers

      DEFAULT_NUMERIC_INCREMENT = '.n'
      SEPARATOR_TOKENS = ['SEP']

      def model_name(identifier)
        matching_rule = rules.find do |rule_name, rule_definition|
          identifier =~ rule_definition.regex
        end

        raise Magma::Gnomon::UnrecognizedIdentifierError.new("#{identifier} does not match any rules.") if matching_rule.nil?

        matching_rule.first
      end

      class << self
        def for_project(project_name)
          self.where(project_name: project_name)
            .reverse(:version_number)
            .first
        end

        def validate(config)
          validation = Magma::Gnomon::Validation.new(config)

          validation.valid? ? [] : validation.errors
        end
      end

      def to_hash
        {
          project_name: project_name,
          config: config,
          version_number: version_number,
          created_at: created_at.iso8601
        }
      end

      def tokens
        config['tokens']
      end

      def rules
        config['rules'].keys.map do |rule_name|
          [ rule_name, rule_parser.fetch(rule_name) ]
        end.to_h
      end

      def decompose(identifier)
        name, rule = rules.find do |name,r| identifier =~ r.regex end

        return nil if !rule

        tokens = rule.decomposition(identifier)

        decomposition = rule_parser.expand_synonyms(tokens.to_h)

        {
          tokens: tokens,
          rules: rules.map do |name, rule|
            rule.from_decomposition(decomposition, project_name)
          end.compact.to_h
        }
      end

      private

      def rule_parser
        @rule_parser ||= Magma::Gnomon::GrammarParser.new(config)
      end
    end
  end
end
