require_relative 'subselect'

class Magma
  class NestingSubselect < Magma::Subselect
    attr_reader :attribute_alias

    def initialize(
      incoming_alias:,
      attribute_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      outgoing_fk_column_name:,
      restrict:,
      subselect:
    )
      super

      @subselect = subselect
    end

    private

    def subselect_data
      identifier_tuple
    end

    def identifier_tuple
      # Returns a tuple like [ ::identifier, <requested value> ]
      Sequel.function(:json_build_array,
        outgoing_identifier_column,
        @subselect
      )
    end
  end
end
