require_relative 'subselect'

class Magma
  class SubselectChild < Magma::Subselect

    private

    def inner_select_contents
      subselect_data
    end
  end
end
