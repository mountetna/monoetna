require_relative 'column_subselect'
require_relative '../with_file_predicate_methods'

class Magma
  class FileSubselectPredicate < Magma::ColumnSubselectPredicate
    include WithFilePredicateMethods

    def self.verbs
      Magma::FilePredicate.verbs
    end

    def select
      []
    end
  end
end
