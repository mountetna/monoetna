require_relative '../match'
require_relative 'with_subselect_override_methods'

class Magma
  class MatchSubselectPredicate < Magma::MatchPredicate
    include WithSubselectOverrideMethods

    def self.verbs
      Magma::MatchPredicate.verbs
    end
  end
end
