require_relative 'column_subselect'

class Magma
  class StringSubselectPredicate < Magma::ColumnSubselectPredicate
    def self.verbs
      Magma::StringPredicate.verbs
    end
  end
end
