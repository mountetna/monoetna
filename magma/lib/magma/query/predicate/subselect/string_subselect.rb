require_relative 'with_subselect_override_methods'
require_relative '../string'

class Magma
  class StringSubselectPredicate < Magma::StringPredicate
    include WithSubselectOverrideMethods

    def self.verbs
      Magma::StringPredicate.verbs
    end
  end
end
