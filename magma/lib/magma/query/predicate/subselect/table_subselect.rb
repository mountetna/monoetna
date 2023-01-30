require_relative 'with_subselect_override_methods'
require_relative '../table'

class Magma
  class TableSubselectPredicate < Magma::TablePredicate
    include WithSubselectOverrideMethods

    def self.verbs
      Magma::TablePredicate.verbs
    end
  end
end
