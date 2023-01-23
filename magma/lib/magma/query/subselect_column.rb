class Magma
  class SubselectColumn
    attr_reader :outgoing_data_column_name, :subselect_column_alias

    def initialize(subselect_column_alias, outgoing_data_column_name)
      @outgoing_data_column_name = outgoing_data_column_name
      @subselect_column_alias = subselect_column_alias
    end
  end
end
