class Magma
  module Gnomon
    class RuleDefinition
      attr_reader :name

      def initialize(grammar, name, definition)
        @grammar = grammar
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

      def regex(with_increment: true)
        /^#{tokenized_definition(with_increment: with_increment)}$/
      end

      def duplicative_tokens?
        seen_tokens = []

        tokens.each do |token|
          next if @grammar.tokens.is_separator?(token) || @grammar.tokens.is_numeric?(token)

          return true if seen_tokens.include?(token)

          seen_tokens << token
        end

        false
      end

      def incrementable?
        !!(raw.strip =~ /\.n$/)
      end

      def illegal_increment_location?
        return false unless raw =~ /\.n/

        !(raw.strip =~ /\.n$/)
      end

      def has_required_tokens?(decomposition)
        (tokens - decomposition.keys).empty?
      end

      def from_decomposition(decomposition, project_name)
        return nil unless has_required_tokens?(decomposition)

        composed_identifier = compose(decomposition)

        identifier_record = Magma::Gnomon::Identifier.where(
          project_name: project_name,
          rule: @name,
          identifier: composed_identifier
        ).first

        magma_model = Magma.instance.get_project(project_name)&.models[ @name.to_sym ]

        magma_record = magma_model&.where(
          magma_model.identity.name => composed_identifier
        )&.first

        [
          @name,
          {
            name: composed_identifier,
            name_created_at: identifier_record&.created_at&.iso8601,
            record_created_at: magma_record&.created_at&.iso8601
          }
        ]
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

            if @grammar.tokens.is_numeric?(token)
              result << "#{@name}_counter"
              next
            end

            rule_name = token[1..-1]

            unless @grammar.rules[rule_name]
              result << token
              next
            end

            raise Magma::Gnomon::RecursiveRuleError.new("Rule \"#{@name}\" may be recursive! Its token \"#{token}\" appears to lead to circular logic.") if expanded_rules.include?(rule_name)

            expanded_rules << rule_name

            result << @grammar.rules[ rule_name ].expanded_definition( expanded_rules )
          end
        end.join(" ")
      end

      def decomposition(identifier)
        tokens.zip(regex.match(identifier).to_a[1..-1])
      end

      def next(identifier_root)
        raise UnincrementableRuleError.new("Rule \"#{name}\" cannot be incremented.") unless incrementable?

        raise UnrecognizedIdentifierError.new("\"#{identifier_root}\" does not match rule \"#{name}\".") unless regex(with_increment: false) =~ identifier_root

        newest_identifier = latest_identifier(identifier_root)

        if newest_identifier.nil?
          next_value = 1
        else
          decomposition = decomposition(newest_identifier.identifier)
          next_value = decomposition.last.last.to_i + 1
        end

        next_value
      end

      def valid?(identifier)
        !!(regex =~ identifier)
      end

      private

      def latest_identifier(identifier_root)
        Magma::Gnomon::Identifier.where(
          rule: name
        ).where { identifier =~ Regexp.new(identifier_root) }.order(:created_at).last
      end

      def all_rules
        @all_rules ||= @config['rules']
      end

      def raw
        @definition
      end

      def tokenized_definition(with_increment: true)
        tokens = expanded_definition.split(" ")

        tokens.pop if !with_increment && @grammar.tokens.is_numeric?(tokens.last)

        tokens.map do |token|
          if @grammar.tokens.is_numeric?(token)
            "(\\d+)"
          else
            "(#{@grammar.tokens.values(token).join("|")})"
          end
        end.join("")
      end

      def self.convert_placeholder_to_name(placeholder)
        placeholder.sub(".", "")
      end
    end
  end
end
