class Magma
  module Gnomon
    class Grammar < Sequel::Model
      def self.to_schema
        {
          "$schema": "http://json-schema.org/draft-07/schema#",
          title: "Gnomon Grammar",
          description: "This grammar defines a set of tokens and rules which can be used to match or create valid identifiers for a Mount Etna project.",
          definitions: [
            Magma::Gnomon::Tokens,
            Magma::Gnomon::Synonyms,
            Magma::Gnomon::Rules,
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

      one_to_many :identifiers

      DEFAULT_NUMERIC_INCREMENT = '.n'
      SEPARATOR_TOKENS = ['SEP']

      def model_name(identifier)
        matching_rule = parser.rules.for(identifier)

        raise Magma::Gnomon::UnrecognizedIdentifierError.new("#{identifier} does not match any rules.") if matching_rule.nil?

        matching_rule.name
      end

      class << self
        def for_project(project_name)
          self.where(project_name: project_name)
            .reverse(:version_number)
            .first
        end

        def validate(config)
          parser = Magma::Gnomon::Grammar::Parser.new(config)

          validation = Magma::Gnomon::Validation.new(parser)

          validation.valid? ? [] : validation.errors
        end
      end

      def to_revision
        to_hash.slice(:config, :created_at, :version_number).merge(comment: comment)
      end

      def to_hash
        {
          project_name: project_name,
          config: config,
          version_number: version_number,
          created_at: created_at.iso8601
        }
      end

      class Parser
        attr_reader :config

        def initialize(config)
          @config = config
        end

        def tokens
          @tokens ||= Magma::Gnomon::Tokens.new(self)
        end

        def rules
          @rules ||= Magma::Gnomon::Rules.new(self)
        end

        def synonyms
          @synonyms ||= Magma::Gnomon::Synonyms.new(self)
        end

        def numeric_increment
          config['numeric_increment'] || Magma::Gnomon::Grammar::DEFAULT_NUMERIC_INCREMENT
        end
      end

      def parser
        @parser ||= Parser.new(config)
      end

      def decompose(identifier)
        rule = parser.rules.for(identifier)

        return nil if !rule

        tokens = rule.decomposition(identifier)

        decomposition = parser.synonyms.expand(tokens.to_h)

        {
          tokens: tokens,
          rule_name: rule.name,
          rules: parser.rules.map do |name, rule|
            rule.from_decomposition(decomposition, project_name)
          end.compact.to_h
        }
      end
    end
  end
end
