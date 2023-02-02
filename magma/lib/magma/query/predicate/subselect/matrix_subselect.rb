require 'set'

class Magma
  class MatrixSubselectPredicate < Magma::MatrixPredicate
    def self.verbs
      Magma::MatrixPredicate.verbs
    end

    def select(incoming_alias_name, incoming_attribute)
      [ Magma::SubselectColumn.new(column_name, attribute_column_name) ]
    end
  end
end
