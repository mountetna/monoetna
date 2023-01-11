require_relative 'record'

class Magma
  class ChildModelRecordPredicate < Magma::RecordPredicate
    def self.verbs
      Magma::RecordPredicate.verbs
    end

    def generate_subselect(incoming_alias_name, incoming_attribute = nil, should_coalesce: false)
      return child_predicate.select unless @is_subselect

      incoming_attribute = valid_attribute(@arguments[0]) if incoming_attribute.nil?

      # If the attribute is a Magma::ChildAttribute, we want to return
      #   a Magma::Subselect here, so the aliases match
      #   up, otherwise the terminal ColumnPredicate won't
      #   have right intermediate Subselect(s).
      [ Magma::Subselect.new(**child_subselect_params(
        incoming_alias_name,
        incoming_attribute
      )) ]
    end

    private

    def child_subselect(incoming_alias_name, incoming_attribute)
      child_predicate.generate_subselect(incoming_alias_name, incoming_attribute, should_coalesce: false).first
    end

    def child_subselect_params(incoming_alias_name, incoming_attribute)
      {
        incoming_alias: incoming_alias_name,
        outgoing_alias: alias_name,
        outgoing_identifier_column_name: @model.identity.column_name,
        outgoing_fk_column_name: outgoing_attribute.column_name,
        outgoing_model: @model,
        restrict: @question.restrict?,
        requested_data: child_subselect(alias_name,
          child_predicate_incoming_attribute),
        filters: [],
        alias_column: false
      }
    end

    def outgoing_attribute
      # Here, for an explicit ChildAttribute, we only need the ParentAttribute
      #   to identify the outgoing attribute.
      @model.attributes.values.select do |attribute|
        attribute.is_a?(Magma::ParentAttribute)
      end.first
    end

    def child_predicate_incoming_attribute
      valid_attribute(@arguments[0])
    end
  end
end
