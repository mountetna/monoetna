
class Magma
  class Grammar < Sequel::Model
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
        ].map(&:to_schema).reduce(&:merge),
	type: "object",
        properties: {
          "tokens": {
            "$ref": "#/definitions/gnomon_tokens"
          },
          "rules": {
            "$ref": "#/definitions/gnomon_rules"
          }
        },
        additionalProperties: false
      }
    end

    class Validation
      attr_reader :errors

      def initialize(config)
        @config = config
        @errors = []
      end

      def valid?
        schema = JSONSchemer.schema(
          JSON.parse(Magma::Grammar.to_schema.to_json)
        )
        schema.valid?(JSON.parse(@config.to_json))
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
    end
  end
end
