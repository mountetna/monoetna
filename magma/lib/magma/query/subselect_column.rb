class Magma
  class SubselectColumn
    attr_reader :outgoing_data_column_name

    def initialize(outgoing_data_column_name)
      @outgoing_data_column_name = outgoing_data_column_name
    end
  end
end
