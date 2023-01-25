require_relative 'record'

class Magma
  class ChildModelPredicate < Magma::RecordPredicate
    def self.verbs
      Magma::RecordPredicate.verbs
    end

    def attribute_select(incoming_alias_name, incoming_attribute = nil)
      [ child_predicate.select ].flatten
    end
  end
end
