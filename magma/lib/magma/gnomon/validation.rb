class Magma
  module Gnomon
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
          JSON.parse(Magma::Gnomon::Grammar.to_schema.to_json)
        )

        schema_errors = schema.validate(JSON.parse(@config.to_json))

        @errors += schema_errors.map do |error|
          JSONSchemer::Errors.pretty(error)
        end
      end

      def validate_rules
        rule_parser = Magma::Gnomon::GrammarParser.new(@config)

        valid_rules = rule_parser.valid?

        @errors += rule_parser.errors unless valid_rules

        valid_rules
      end

      def validate_synonyms
        # Are synonym labels unique
        # Do all synonyms have equal numbers of values
        # Do all synonyms have matching values
        return true unless @config['synonyms'] && @config['synonyms'].length > 0

        validator = Magma::Gnomon::SynonymValidator.new(@config)

        valid = validator.valid?

        @errors += validator.errors unless valid

        valid
      end
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
          :validate_no_duplicate_tokens,
          :validate_matching_values
        ]
      end

      def validate_no_duplicate_tokens
        counts = @config['synonyms'].flatten.tally.select do |tok, count|
          count > 1
        end

        @errors << "Duplicate tokens #{counts.keys.join(", ")} in synonyms" unless counts.empty?
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
  end
end
