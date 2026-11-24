require_relative 'query_validation'

class Magma
  class SearchField < Sequel::Model
    include Magma::Search::QueryValidation
    plugin :auto_validations, not_null: :presence

    many_to_one :modality, class: 'Magma::SearchModality'

    def to_hash
      { id: id, label: label, query: query.to_a, kind: kind, modality_id: modality_id,
        modality: modality&.name, position: position, updated_at: updated_at&.iso8601 }
    end

    def validate
      super
      check_query!(query)
      return unless (verbs = VERBS[kind]) # a bad kind is left to the table's CHECK

      verb = result_verb(query.to_a)
      return if verbs.include?(verb) && (kind == 'count' || query.to_a.last != verb)

      raise Etna::BadRequest, 'A count must end in ::count. A filter must end in ::distinct or ::first, followed by an attribute.'
    end

    def before_create
      self.position ||= (self.class.max(:position) || 0) + 1
      super
    end
  end
end
