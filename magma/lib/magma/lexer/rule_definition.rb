class Magma
  class RuleDefinition
    attr_reader :name

    def initialize(parser, name, definition)
      @parser = parser
      @name = name
      @definition = definition
    end

    def empty?
      raw.strip.empty?
    end

    def placeholder(name=nil)
      ".#{name || @name}"
    end

    def tokens
      expanded_definition.split(" ")
    end

    def regex
      /^#{tokenized_definition}$/
    end

    def duplicative_tokens?
      seen_tokens = []

      tokens.each do |token|
        next if @parser.is_separator_token?(token) || @parser.is_numeric_token?(token)

        return true if seen_tokens.include?(token)

        seen_tokens << token
      end

      false
    end

    def illegal_increment_location?
      return false unless raw =~ /\.n/

      !(raw.strip =~ /\.n$/)
    end

    def has_required_tokens?(decomposition)
      (tokens - decomposition.keys).empty?
    end

    def compose(decomposition)
      tokens.map do |token|
        decomposition[token]
      end.join('')
    end

    def expanded_definition(expanded_rules = nil)
      @expanded_definition ||= [].tap do |result|
        expanded_rules ||= Set.new([ @name ])

        raw.split(" ").each do |token|
          if token[0] != '.'
            result << token
            next
          end

          if @parser.is_numeric_token?(token)
            result << "#{@name}_counter"
            next
          end

          rule_name = token[1..-1]

          unless @parser.rules[rule_name]
            result << token
            next
          end

          raise Magma::RecursiveRuleError.new("Rule \"#{@name}\" may be recursive! Its token \"#{token}\" appears to lead to circular logic.") if expanded_rules.include?(rule_name)

          expanded_rules << rule_name

          result << @parser.rules[ rule_name ].expanded_definition( expanded_rules )
        end
      end.join(" ")
    end

    def decomposition(identifier)
      tokens.zip(regex.match(identifier).to_a[1..-1])
    end

    private

    def all_rules
      @all_rules ||= @config['rules']
    end

    def all_tokens
      @parser.tokens
    end

    def raw
      @definition
    end

    def tokenized_definition
      expanded_definition.split(" ").map do |token|
        if @parser.is_numeric_token?(token)
          "(\\d+)"
        else
          "(#{values_for_token(token).join("|")})"
        end
      end.join("")
    end

    def values_for_token(token)
      all_tokens[token]&.fetch('values')&.keys || []
    end

    def self.convert_placeholder_to_name(placeholder)
      placeholder.sub(".", "")
    end
  end
end
