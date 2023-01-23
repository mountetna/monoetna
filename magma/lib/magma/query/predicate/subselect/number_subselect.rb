require_relative 'column_subselect'

class Magma
  class NumberSubselectPredicate < Magma::ColumnSubselectPredicate
    def self.verbs
      Magma::NumberPredicate.verbs
    end
  end
end
