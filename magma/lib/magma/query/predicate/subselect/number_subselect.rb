require_relative 'with_subselect_override_methods'
require_relative '../number'

class Magma
  class NumberSubselectPredicate < Magma::NumberPredicate
    include WithSubselectOverrideMethods

    def self.verbs
      Magma::NumberPredicate.verbs
    end
  end
end
