class Magma
  module Gnomon
    class Synonyms
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

      attr_reader :errors

      def initialize(grammar)
        @grammar = grammar
        @errors = []
      end

      def config
        @grammar.config['synonyms'] || []
      end

      def valid?
        return true unless config && !config.empty?

        validations.each do |validation|
          send(validation)
        end

        @errors.empty?
      end

      def expand(decomposition)
        # add any synonym tokens to decomposition
        expansion = {}
        decomposition.each do |token_name, token_value|
          next if @grammar.tokens.is_separator?(token_name)

          config.each do |syn_set|
            next unless syn_set.include?(token_name)
            syn_set.each do |synonym|
              expansion[synonym] = @grammar.tokens.match_value_by_description(synonym, token_name, token_value)
            end
          end
        end
        decomposition.merge(expansion)
      end

      private

      def validations
        [
          :validate_no_undefined_tokens,
          :validate_unique_values,
          :validate_no_duplicate_tokens,
          :validate_matching_values
        ]
      end

      def validate_no_duplicate_tokens
        counts = duplicates(config.flatten)

        @errors << "Duplicate tokens #{counts.join(", ")} in synonyms" unless counts.empty?
      end

      def validate_no_undefined_tokens
        config.each do |synonym_set|
          missing_tokens = synonym_set - @grammar.tokens.valid_token_names(synonym_set)

          @errors << "Missing token definitions for: #{comma_separate_list(missing_tokens)}, for synonyms #{comma_separate_list(synonym_set)}." if missing_tokens.length > 0
        end
      end

      def validate_unique_values
        config.each do |synonym_set|
          tokens = @grammar.tokens.valid_tokens(synonym_set)

          nonunique_values = tokens.any? do |token|
            !duplicates(token['values'].values).empty?
          end

          @errors << "Synonyms #{comma_separate_list(synonym_set)} do not have unique values." if nonunique_values
        end
      end

      def validate_matching_values
        config.each do |synonym_set|
          token_descs = @grammar.tokens.valid_descriptions(synonym_set)

          next if token_descs.empty?

          ref = token_descs.first

          different_values = token_descs.any? do |descs|
            !((ref-descs) + (descs-ref)).empty?
          end

          @errors << "Synonyms #{comma_separate_list(synonym_set)} do not have matching values." if different_values
        end
      end

      def duplicates(array)
        array.tally.select do |o, count|
          count > 1
        end.keys
      end

      def comma_separate_list(strings)
        strings.sort.map { |t| "\"#{t}\""}.join(", ")
      end
    end
  end
end
