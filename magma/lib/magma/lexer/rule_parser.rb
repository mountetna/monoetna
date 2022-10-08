class Magma
  class RuleParser
    attr_reader :errors

    def initialize(config)
      @config = config
      @errors = []
    end

    def valid?
      validations.each do |validation|
        send(validation)
      end

      @errors.empty?
    end

    def fetch(rule_name)
      raise UndefinedRuleError.new("#{rule_name} undefined") unless rules[rule_name]

      rules[rule_name]
    end

    def rules
      @rules ||= @config['rules']&.map do |rule_name, rule_definition|
        [
          rule_name, 
          Magma::RuleDefinition.new(
            self,
            rule_name,
            rule_definition
          )
        ]
      end.to_h || {}
    end

    def rule_names
      rules.keys
    end

    def synonyms
      @synonyms ||= @config['synonyms'] || []
    end

    def tokens
      @tokens ||= @config['tokens'] || {}
    end

    def expand_synonyms(decomposition)
      # add any synonym tokens to decomposition

      merge = {}
      decomposition.each do |token_name, token_value|
        next if Magma::Grammar::SEPARATOR_TOKENS.include?(token_name)

        synonyms.each do |syn_set|
          next unless syn_set.include?(token_name)
          syn_set.each do |syn|
            key, value = tokens[syn]['values'].find do |key,value|
              value == tokens[token_name]['values'][token_value]
            end
            merge[syn] = key
          end
        end
      end
      decomposition.merge(merge)
    end

    def is_separator_token?(token)
      Magma::Grammar::SEPARATOR_TOKENS.include?(token)
    end

    def is_numeric_token?(token)
      (@config['numeric_increment'] || Magma::Grammar::DEFAULT_NUMERIC_INCREMENT) == token || token =~ /_counter$/
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
        :validate_single_separator
      ]
    end

    def numeric_increment
      @config['numeric_increment'] || Magma::Grammar::DEFAULT_NUMERIC_INCREMENT
    end

    def validate_no_duplicate_rules
      seen_rule_definitions = []
      duplicate_definitions = []

      # check using the regex to avoid ambiguity
      expanded_rules = rules.map do |rule_name, rule_definition|
        [ rule_name, rule_definition.regex ]
      rescue Magma::RecursiveRuleError => e
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

    def all_token_placeholders
      @all_token_values ||= (@config['tokens']&.keys || [])
    end

    def all_rule_placeholders
      @all_rule_placeholders ||= rules.map do |rule_name, rule_definition|
        rule_definition.placeholder
      end
    end

    def allowed_rule_tokens
      @allowed_rule_tokens ||= all_token_placeholders + all_rule_placeholders + [numeric_increment]
    end

    def validate_all_tokens_defined
      unidentified_tokens = []

      unknown_token_exists = @config['rules']&.any? do |rule_name, rule_definition|
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
      rescue Magma::RecursiveRuleError => e
        recursive_errors << e.message
      end

      @errors += recursive_errors unless recursive_errors.empty?

      recursive_errors.empty?
    end

    def validate_no_duplicate_tokens
      duplicate_tokens_errors = []

      rules.each do |rule_name, rule_definition|
        duplicate_tokens_errors << "Rule \"#{rule_name}\" contains duplicate tokens." if rule_definition.duplicative_tokens?
      rescue Magma::RecursiveRuleError => e
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

    def validate_single_separator
      all_tokens = @config['tokens'] || {}

      defined_separators = Magma::Grammar::SEPARATOR_TOKENS & all_tokens.keys
      @errors << "No separator token defined!" unless defined_separators.size > 0
      @errors << "More than one separator token defined!" unless defined_separators.size < 2
      @errors << "More than one separator token defined!" unless defined_separators.all? {|sep| all_tokens[sep].fetch('values').size == 1}
    end
  end
end
