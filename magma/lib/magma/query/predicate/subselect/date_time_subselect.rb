require_relative '../date_time'
require_relative 'with_subselect_override_methods'

class Magma
  class DateTimeSubselectPredicate < Magma::DateTimePredicate
    include WithSubselectOverrideMethods

    def self.verbs
      Magma::DateTimePredicate.verbs
    end
  end
end
