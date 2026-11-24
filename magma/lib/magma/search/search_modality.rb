require_relative 'query_validation'

class Magma
  class SearchModality < Sequel::Model
    include Magma::Search::QueryValidation
    plugin :auto_validations, not_null: :presence

    one_to_many :fields, class: 'Magma::SearchField', key: :modality_id

    def to_hash
      { id: id, name: name, model: self[:model], position: position, updated_at: updated_at&.iso8601 }
    end

    # The unique index is on lower(name); trimming keeps ' Flow' and 'flow' the same.
    def before_validation
      self.name = name.to_s.strip
      super
    end

    def validate
      super
      check_template_model!(self[:model])
    end

    def before_create
      self.position ||= (self.class.max(:position) || 0) + 1
      super
    end

    def before_destroy
      raise Etna::BadRequest, "'#{name}' is used by #{fields_dataset.count} field(s); move or delete them first" if fields_dataset.any?
      super
    end
  end
end
