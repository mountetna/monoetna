require_relative 'with_table_predicate_methods'

class Magma
  # just like the ModelPredicate, this keeps track of its own predicate chain.
  # Confusing... Perhaps a better concept is in order?
  class TablePredicate < Magma::Predicate
    include WithTablePredicateMethods

    verb nil do
      child Array
    end
  end
end
