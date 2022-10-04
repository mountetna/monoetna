require_relative './lexer/rule_parser'

class Magma
  class Grammar < Sequel::Model
    class Synonym
      def self.to_schema
        {
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
      end
    end
    class Rule
      def self.to_schema
        {
          gnomon_rules: {
            type: "object",
            patternProperties: {
              "^\\w+$": {
                type: "string"
              }
            }
          }
        }
      end
    end
    class Token
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
                }
              }
            }
          }
        }
      end
    end
    def self.to_schema
      {
        "$schema": "http://json-schema.org/draft-07/schema#",
	title: "Gnomon Grammar",
        description: "This grammar defines a set of tokens and rules which can be used to match or create valid identifiers for a Mount Etna project.",
        definitions: [
          Magma::Grammar::Rule,
          Magma::Grammar::Token,
          Magma::Grammar::Synonym,
        ].map(&:to_schema).reduce(&:merge),
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

    class Validation
      def initialize(config)
        @config = config
      end

      def valid?
        schema = JSONSchemer.schema(
          JSON.parse(Magma::Grammar.to_schema.to_json)
        )

        errors = schema.validate(JSON.parse(@config.to_json))

        errors.map do |error|
          JSONSchemer::Errors.pretty(error)
        end
      end
    end

    one_to_many :identifiers

    class << self
      def for_project(project_name)
        self.where(project_name: project_name)
          .reverse(:version_number)
          .first
      end

      def validate(config)
        Magma::Grammar::Validation.new(config).valid?
      end
    end

    def to_hash
      {
        project_name: project_name,
        config: config,
        version_number: version_number,
        created_at: created_at
      }
    def tokens
      config['tokens']
    end

    def rules
      config['rules'].map do |token, rule|
        rule_parser.parse(rule)
        # Then what happens??
        Magma::Rule.new(token: token, regex: rule_parser.parse(rule))
      end
    end

    private

    def rule_parser
      @rule_parser || Magma::RuleParser.new(self)
    end
  end
end
