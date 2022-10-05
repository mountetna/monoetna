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

    class GrammarError < Exception
    end

    class UnrecognizedIdentifierError < GrammarError
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

    class SynonymValidator
      attr_reader :errors

      def initialize(config)
        @config = config
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
          :validate_no_undefined_tokens,
          :validate_tokens_same_length,
          :validate_unique_values,
          :validate_matching_values
        ]
      end

      def validate_no_undefined_tokens
        @config['synonyms'].each do |synonym_set|
          missing_tokens = synonym_set - tokens_for_synonyms(synonym_set).keys

          @errors << "Missing token definitions for: #{comma_separate_list(missing_tokens)}, for synonyms #{comma_separate_list(synonym_set)}." if missing_tokens.length > 0
        end
      end

      def tokens_for_synonyms(synonym_set)
        @config['tokens'].dup.slice(*synonym_set)
      end

      def validate_tokens_same_length
        @config['synonyms'].each do |synonym_set|
          tokens = tokens_for_synonyms(synonym_set)
          reference_length = tokens.values.first['values'].keys.length
          different_lengths = tokens.any? do |token_name, token_definition|
            token_definition['values'].keys.length != reference_length
          end

          @errors << "Synonyms #{comma_separate_list(synonym_set)} do not have an equal number of values." if different_lengths
        end
      end

      def validate_unique_values
        @config['synonyms'].each do |synonym_set|
          tokens = tokens_for_synonyms(synonym_set)
          nonunique_values = tokens.any? do |token_name, token_definition|
            seen_values = []
            token_definition['values'].values.any? do |value|
              seen_values.include?(value) ? true : begin
                seen_values << value
                false
              end
            end
          end

          @errors << "Synonyms #{comma_separate_list(synonym_set)} do not have unique values." if nonunique_values
        end
      end

      def validate_matching_values
        @config['synonyms'].each do |synonym_set|
          tokens = tokens_for_synonyms(synonym_set)
          reference_values = Set.new(tokens.values.first['values'].values)
          different_values = tokens.any? do |token_name, token_definition|
            Set.new(token_definition['values'].values) != reference_values
          end

          @errors << "Synonyms #{comma_separate_list(synonym_set)} do not have matching values." if different_values
        end
      end

      def comma_separate_list(strings)
        strings.sort.map { |t| "\"#{t}\""}.join(", ")
      end
    end

    class Validation
      attr_reader :errors

      def initialize(config)
        @config = config
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
          :validate_synonyms
        ]
      end

      def validate_schema
        schema = JSONSchemer.schema(
          JSON.parse(Magma::Grammar.to_schema.to_json)
        )

        schema_errors = schema.validate(JSON.parse(@config.to_json))

        @errors += schema_errors.map do |error|
          JSONSchemer::Errors.pretty(error)
        end
      end

      def validate_rules
        rule_parser = Magma::RuleParser.new(@config)

        valid_rules = rule_parser.valid?

        @errors += rule_parser.errors unless valid_rules

        valid_rules
      end

      def validate_synonyms
        # Are synonym labels unique
        # Do all synonyms have equal numbers of values
        # Do all synonyms have matching values
        return true unless @config['synonyms'] && @config['synonyms'].length > 0

        validator = Magma::Grammar::SynonymValidator.new(@config)

        valid = validator.valid?

        @errors += validator.errors unless valid

        valid
      end
    end

    one_to_many :identifiers

    DEFAULT_NUMERIC_INCREMENT = '.n'
    SEPARATOR_TOKENS = ['SEP', 'SEPARATOR']

    def tokens
      config['tokens'].map do |token_name, token_definition|
        Magma::Grammar::Token.new(token_name, token_definition)
      end
    end

    def rules
      @rules ||= config['rules'].keys.map do |rule_name|
        [rule_name, rule_parser.fetch(rule_name)]
      end.to_h
    end

    def model_name(identifier)
      matching_rule = rules.find do |rule_name, rule_definition|
        identifier =~ rule_definition.regex
      end

      raise UnrecognizedIdentifierError.new("#{identifier} does not match any rules.") if matching_rule.nil?

      matching_rule.first
    end

    class << self
      def for_project(project_name)
        self.where(project_name: project_name)
          .reverse(:version_number)
          .first
      end

      def validate(config)
        validation = Magma::Grammar::Validation.new(config)

        validation.valid? ? [] : validation.errors
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

    def tokens
      config['tokens']
    end

    def rules
      config['rules'].keys.map do |rule_name|
        [ rule_name, rule_parser.fetch(rule_name) ]
      end.to_h
    end

    private

    def rule_parser
      @rule_parser ||= begin
        parser = Magma::RuleParser.new(config)

        valid = parser.valid?
        raise Magma::RuleParser::RuleParseError.new("Config errors: #{parser.errors}") unless valid

        parser
      end
    end
  end
end
