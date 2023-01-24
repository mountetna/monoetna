require_relative 'column_subselect'
require_relative '../with_date_time_predicate_methods'

class Magma
  class DateTimeSubselectPredicate < Magma::ColumnSubselectPredicate
    include WithDateTimePredicateMethods

    def self.verbs
      Magma::DateTimePredicate.verbs
    end
  end
end
