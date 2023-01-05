require_relative 'subselect'

class Magma
  class TerminalSubselect < Magma::Subselect
    private

    def subselect_data
      outgoing_column
    end
  end
end
