require_relative 'with_date_time_predicate_methods'

class Magma
  class DateTimePredicate < Magma::ColumnPredicate
    include WithDateTimePredicateMethods

    verb nil do
      child DateTime
    end

    verb [ '::<=', '::<', '::>=', '::>', '::=', '::!=' ], String do
      child TrueClass

      constraint do
        op, date = @arguments
        comparison_constraint(@column_name, op, DateTime.parse(date))
      end
    end
  end
end
