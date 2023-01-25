require_relative 'column_subselect'
require_relative '../with_table_predicate_methods'

class Magma
  class TableSubselectPredicate < Magma::ColumnSubselectPredicate
    include WithTablePredicateMethods

    def self.verbs
      Magma::TablePredicate.verbs
    end
  end
end
