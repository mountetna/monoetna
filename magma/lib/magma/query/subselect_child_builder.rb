require_relative 'subselect_builder'

class Magma
  class SubselectChildBuilder < Magma::SubselectBuilder

    private

    def inner_select_contents
      subselect_data
    end
  end
end
