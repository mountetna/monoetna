require_relative '../boolean'
require_relative 'with_subselect_override_methods'

class Magma
  class BooleanSubselectPredicate < Magma::BooleanPredicate
    include WithSubselectOverrideMethods

    def self.verbs
      Magma::BooleanPredicate.verbs
    end
  end
end
