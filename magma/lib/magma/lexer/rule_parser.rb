class Magma
  class RuleParser
    class RuleParseError < Exception
    end

    class RecursiveRuleError < RuleParseError
    end

    class UndefinedRuleError < RuleParseError
    end

    class RuleDefinition
      attr_reader :name

      def initialize(name, config, expansion_cache={})
        @name = name
        @config = config
        @expansion_cache = expansion_cache
      end

      def self.from_placeholder(placeholder, config, expansion_cache)
        rule_name = self.convert_placeholder_to_name(placeholder)

        Magma::RuleParser::RuleDefinition.new(
          rule_name,
          config,
          expansion_cache
        )
      end

      def empty?
        raw.strip.empty?
      end

      def placeholder
        ".#{@name}"
      end

      def regex(with_increment: true)
        /^#{tokenized_definition(with_increment: with_increment)}$/
      end

      def incrementable?
        raw.strip =~ /\.n$/
      end

      def duplicative_tokens?
        seen_tokens = []

        expanded_definition.split(" ").each do |token|
          next if Magma::Grammar::SEPARATOR_TOKENS.include?(token) || numeric_increment == token

          return true if seen_tokens.include?(token)

          seen_tokens << token
        end

        false
      end

      def illegal_increment_location?
        return false unless raw =~ /\.n/

        !(raw.strip =~ /\.n$/)
      end

      def expanded_definition(seen_placeholders = [])
        @expanded_definition ||= [].tap do |result|
          seen_placeholders << placeholder

          raw.split(" ").map do |token|

            raise RecursiveRuleError.new("Rule \"#{@name}\" may be recursive! It's token \"#{token}\" appears to lead to circular logic.") if seen_placeholders.include?(token)

            seen_placeholders << token if other_rule_placeholders.include?(token)

            if other_rule_placeholders.include?(token) && @expansion_cache.include?(token)
              result << @expansion_cache[token]
            elsif other_rule_placeholders.include?(token)
              other_definition = RuleDefinition.from_placeholder(
                token,
                @config,
                @expansion_cache
              )

              other_expanded_definition = other_definition.expanded_definition(seen_placeholders)
              @expansion_cache[token] = other_expanded_definition
              result << other_expanded_definition
            else
              result << token
            end
          end
        end.join(" ")
      end

      private

      def all_rules
        @all_rules ||= @config['rules']
      end

      def all_tokens
        @all_tokens ||= @config['tokens']
      end

      def raw
        @raw ||= all_rules[@name]
      end

      def tokenized_definition(with_increment: true)
        tokens = expanded_definition.split(" ")

        tokens.pop if !with_increment && tokens.last == numeric_increment

        tokens.map do |token|
          if token == numeric_increment
            "\\d+"
          else
            "(#{values_for_token(token).join("|")})"
          end
        end.join("")
      end

      def values_for_token(token)
        all_tokens[token]['values'].keys
      end

      def self.convert_placeholder_to_name(placeholder)
        placeholder.sub(".", "")
      end

      def other_rule_placeholders
        other_names = all_rules.keys.reject do |rule_name|
          rule_name == @name
        end.map do |other_name|
          ".#{other_name}"
        end
      end

      def numeric_increment
        @config['numeric_increment'] || Magma::Grammar::DEFAULT_NUMERIC_INCREMENT
      end
    end

    attr_reader :errors

    def initialize(config)
      @config = config
      @rules = {}
      @errors = []
    end

    def valid?
      construct_rules

      validations.each do |validation|
        send(validation)
      end

      @errors.empty?
    end

    def fetch(rule_name)
      raise UndefinedRuleError.new("#{rule_name} undefined") unless @rules[rule_name]

      @rules[rule_name]
    end

    def construct_rules
      expansion_cache = {}
      @config['rules']&.map do |rule_name, rule_definition|
        @rules[rule_name] = RuleDefinition.new(
          rule_name,
          @config,
          expansion_cache
        )
      end
    end

    private

    def validations
      [
        :validate_no_duplicate_rules,
        :validate_no_blank_rules,
        :validate_all_tokens_defined,
        :validate_no_recursive_rules,
        :validate_no_duplicate_tokens,
        :validate_numeric_increment_at_end,
      ]
    end

    def numeric_increment
      @config['numeric_increment'] || Magma::Grammar::DEFAULT_NUMERIC_INCREMENT
    end

    def validate_no_duplicate_rules
      seen_rule_definitions = []
      duplicate_definitions = []

      # First get RuleDefinition for all rules, and
      #   collect the expanded_definition. Then
      #   check for duplicates.
      expanded_rules = @rules.map do |rule_name, rule_definition|
        rule_definition.expanded_definition
      rescue RecursiveRuleError => e
        nil # We'll validate recursive rules in a different validation
      end.compact

      duplicate_exists = expanded_rules.any? do |rule_definition|
        duplicated = seen_rule_definitions.include?(rule_definition)

        (duplicated ? duplicate_definitions : seen_rule_definitions) << rule_definition

        duplicated
      end

      @errors << "Duplicate rule definition exists: #{duplicate_definitions}" if duplicate_exists

      !duplicate_exists
    end

    def validate_no_blank_rules
      blank_rule = @rules.any? do |rule_name, rule_definition|
        rule_definition.empty?
      end

      @errors << "Rules cannot contain only whitespace" if blank_rule

      blank_rule
    end

    def all_token_placeholders
      @all_token_values ||= (@config['tokens']&.keys || [])
    end

    def all_rule_placeholders
      @all_rule_placeholders ||= @rules.map do |rule_name, rule_definition|
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

      @rules.each do |rule_name, rule_definition|
        rule_definition.expanded_definition
      rescue RecursiveRuleError => e
        recursive_errors << e.message
      end

      @errors += recursive_errors unless recursive_errors.empty?

      recursive_errors.empty?
    end

    def validate_no_duplicate_tokens
      duplicate_tokens_errors = []

      @rules.each do |rule_name, rule_definition|
        duplicate_tokens_errors << "Rule \"#{rule_name}\" contains duplicate tokens." if rule_definition.duplicative_tokens?
      rescue RecursiveRuleError => e
        nil # We'll validate recursive rules in a different validation
      end

      @errors += duplicate_tokens_errors unless duplicate_tokens_errors.empty?

      duplicate_tokens_errors.empty?
    end

    def validate_numeric_increment_at_end
      numeric_increment_at_end_errors = []

      @rules.each do |rule_name, rule_definition|
        numeric_increment_at_end_errors << "Rule \"#{rule_name}\" can only use the numeric increment \".n\" at the end." if rule_definition.illegal_increment_location?
      end

      @errors += numeric_increment_at_end_errors unless numeric_increment_at_end_errors.empty?

      numeric_increment_at_end_errors.empty?
    end
  end
end
