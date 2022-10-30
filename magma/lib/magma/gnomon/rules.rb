class Magma
  module Gnomon
    class Rules
      include Enumerable
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

      def fetch(rule_name)
        raise Magma::Gnomon::UndefinedRuleError.new("#{rule_name} undefined") unless rules[rule_name]

        rules[rule_name]
      end

      def for(identifier)
        name, rule = rules.find do |rule_name, rule_definition|
          identifier =~ rule_definition.regex
        end

        return rule
      end

      def each
        rules.each do |name, rule|
          yield name, rule
        end
      end

      def [](rule_name)
        rules[rule_name]
      end

      def rule_names
        rules.keys
      end

      private

      def validations
        [
          :validate_all_tokens_defined,
          :validate_no_duplicate_rules,
          :validate_no_blank_rules,
          :validate_no_recursive_rules,
          :validate_no_duplicate_tokens,
          :validate_numeric_increment_at_end,
        ]
      end

      def validate_no_duplicate_rules
        seen_rule_definitions = []
        duplicate_definitions = []

        # check using the regex to avoid ambiguity
        expanded_rules = rules.map do |rule_name, rule_definition|
          [ rule_name, rule_definition.regex ]
        rescue Magma::Gnomon::RecursiveRuleError => e
          nil # We'll validate recursive rules in a different validation
        end.compact.to_h

        duplicate_exists = expanded_rules.keys.group_by do |rule_name|
          expanded_rules[rule_name]
        end.select do |rule_def, rule_names|
          rule_names.count > 1
        end

        @errors << "Duplicate rule definition exists: #{duplicate_exists.map {|rule_def,names| "[ #{names.join(', ')} ]"}.join(", ")}" unless duplicate_exists.empty?

        duplicate_exists.empty?
      end

      def validate_no_blank_rules
        blank_rule = rules.any? do |rule_name, rule_definition|
          rule_definition.empty?
        end

        @errors << "Rules cannot contain only whitespace" if blank_rule

        blank_rule
      end

      def allowed_rule_tokens
        @allowed_rule_tokens ||= begin
          all_token_placeholders = @grammar.tokens.token_names

          all_rule_placeholders = rules.map do |rule_name, rule_definition|
            rule_definition.placeholder
          end

          all_token_placeholders + all_rule_placeholders + [@grammar.numeric_increment]
        end
      end

      def validate_all_tokens_defined
        unidentified_tokens = []

        unknown_token_exists = config&.any? do |rule_name, rule_definition|
          rule_definition.split(" ").any? do |rule_token|
            unknown = !allowed_rule_tokens.include?(rule_token)

            unidentified_tokens << rule_token if unknown

            unknown
          end
        end

        @errors << "Unknown tokens used: #{unidentified_tokens}" if unknown_token_exists

        !unknown_token_exists
      end

      def validate_no_recursive_rules
        recursive_errors = []

        rules.each do |rule_name, rule_definition|
          rule_definition.expanded_definition
        rescue Magma::Gnomon::RecursiveRuleError => e
          recursive_errors << e.message
        end

        @errors += recursive_errors unless recursive_errors.empty?

        recursive_errors.empty?
      end

      def validate_no_duplicate_tokens
        duplicate_tokens_errors = []

        rules.each do |rule_name, rule_definition|
          duplicate_tokens_errors << "Rule \"#{rule_name}\" contains duplicate tokens." if rule_definition.duplicative_tokens?
        rescue Magma::Gnomon::RecursiveRuleError => e
          nil # We'll validate recursive rules in a different validation
        end

        @errors += duplicate_tokens_errors unless duplicate_tokens_errors.empty?

        duplicate_tokens_errors.empty?
      end

      def validate_numeric_increment_at_end
        numeric_increment_at_end_errors = []

        rules.each do |rule_name, rule_definition|
          numeric_increment_at_end_errors << "Rule \"#{rule_name}\" can only use the numeric increment \".n\" at the end." if rule_definition.illegal_increment_location?
        end

        @errors += numeric_increment_at_end_errors unless numeric_increment_at_end_errors.empty?

        numeric_increment_at_end_errors.empty?
      end

      def config
        @grammar.config['rules'] || {}
      end

      def rules
        @rules ||= config&.map do |rule_name, rule_definition|
          [
            rule_name, 
            Magma::Gnomon::RuleDefinition.new(
              @grammar,
              rule_name,
              rule_definition
            )
          ]
        end.to_h || {}
      end
    end
  end
end
