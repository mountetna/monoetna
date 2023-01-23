require_relative 'column_subselect'

class Magma
  class MatchSubselectPredicate < Magma::ColumnSubselectPredicate
    def self.verbs
      Magma::MatchPredicate.verbs
    end
  end
end
