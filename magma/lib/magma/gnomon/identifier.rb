class Magma
  module Gnomon
    class Identifier < Sequel::Model
      many_to_one :grammar
      many_to_one :renamed_to, class: self, key: :renamed_to_id
      one_to_many :aliases, class: self, key: :renamed_to_id

      def to_hash(user)
        {
          identifier: identifier,
          name_created_at: created_at,
          author: author,
          record_created_at: record_created_at(user)
        }
      end

      private

      def record_created_at(user)
        record(user)&.first&.last
      end

      def record(user)
        model = Magma.instance.get_model(project_name, rule)

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
    end
  end
end
