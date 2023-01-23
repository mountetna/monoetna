require_relative 'column_subselect'

class Magma
  class BooleanSubselectPredicate < Magma::ColumnSubselectPredicate
    def self.verbs
      Magma::BooleanPredicate.verbs
    end
  end
end
