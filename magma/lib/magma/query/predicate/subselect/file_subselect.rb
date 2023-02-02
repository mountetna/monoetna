require_relative 'with_subselect_override_methods'
require_relative '../file'

class Magma
  class FileSubselectPredicate < Magma::FilePredicate
    include WithSubselectOverrideMethods

    def self.verbs
      Magma::FilePredicate.verbs
    end

    def select
      []
    end
  end
end
