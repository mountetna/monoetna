class Magma
  module Gnomon
    class Identifier < Sequel::Model
      many_to_one :grammar
      many_to_one :renamed_to, class: self, key: :renamed_to_id
      one_to_many :aliases, class: self, key: :renamed_to_id

      def to_hash
        {
          identifier: identifier,
          name_created_at: created_at.iso8601,
          author: author
        }
      end

      private

      def record_created_at(user)
        record(user)&.first&.last
      end

      def record(user)
        question = Magma::Question.new(
          project_name,
          record_query,
          show_disconnected: false,
          restrict: !user.can_see_restricted?(project_name),
          user: user,
          timeout: Magma.instance.config(:query_timeout)
        )

        question.answer
      end

      def record_query
        [rule, ["::identifier", "::equals", identifier], "::all", "created_at"]
      end

      def self.backfill(project_name, rule_name, author, dry_run: true)
        grammar = Magma::Gnomon::Grammar.for_project(project_name)
        raise "No grammar for #{project_name}" unless grammar
        rule = grammar.parser.rules[rule_name]
        raise "No such #{rule_name} for #{project_name}" unless rule
        
        magma_model = Magma.instance.get_model(project_name, rule_name)
        identifiers = magma_model.select_map(magma_model.identity.attribute_name.to_sym)
        bad_identifiers = []
        good_identifiers = []
        identifiers.each do |identifier|
          unless rule.valid?(identifier)
            bad_identifiers << identifier
            next
          end
          good_identifiers << {
            project_name: project_name,
            rule: rule_name,
            author: author,
            created_at: DateTime.now,
            identifier: identifier,
            grammar_id: grammar.id
          }
        end

        # Database insertion
        Magma::Gnomon::Identifier.multi_insert(good_identifiers) unless dry_run

        puts "#{good_identifiers.length} ids backfilled, #{bad_identifiers.length} do not match current grammar."
        unless bad_identifiers.empty?
          puts "Bad Identifiers:" 
          puts bad_identifiers.join(", ")
        end
      end
    end
  end
end
