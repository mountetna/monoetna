class Magma
  class SubselectError < Exception
  end

  class Subselect
    attr_reader :subselect_column_alias

    def initialize(
      incoming_alias:,
      subselect_column_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      restrict:,
      outgoing_fk_column_name:,
      requested_data:
    )
      @incoming_alias = incoming_alias.to_sym
      @subselect_column_alias = subselect_column_alias.to_sym
      @outgoing_model = outgoing_model
      @outgoing_alias = outgoing_alias.to_sym
      @outgoing_identifier_column_name = outgoing_identifier_column_name.to_sym
      @restrict = restrict
      @outgoing_fk_column_name = outgoing_fk_column_name
      @requested_data = requested_data # another subselect, or a column mame on the outgoing model
    end

    def coalesce
      Sequel.function(
        :coalesce,
        inner_select,
        default_value
      ).as(subselect_column_alias)
    end

    def hash
      subselect_column_alias.hash
    end

    def eql?(other)
      subselect_column_alias == other.subselect_column_alias
    end

    private

    def inner_select
      subselect_query = @outgoing_model.from(
        Sequel.as(@outgoing_model.table_name, @outgoing_alias)
      )

      subselect_query.select(
        Sequel.function(
          :json_agg,
          subselect_data
        )
      ).where(**subselect_constraints).where(restrict_constraints)
    end

    def subselect_data
      identifier_tuple
    end

    def identifier_tuple
      # Returns a tuple like [ ::identifier, <requested data> ]
      Sequel.function(:json_build_array,
        outgoing_identifier_column,
        outgoing_data
      )
    end

    def outgoing_identifier_column
      Sequel.qualify(@outgoing_alias, @outgoing_identifier_column_name)
    end

    def outgoing_data
      @requested_data.is_a?(Magma::Subselect) ?
        @requested_data :
        Sequel.qualify(@outgoing_alias, @requested_data.outgoing_data_column_name)
    end

    def default_value
      Sequel.cast('[]', :json)
    end

    def subselect_constraints
      {
        Sequel.qualify(@outgoing_alias, @outgoing_fk_column_name) => Sequel.qualify(@incoming_alias, :id)
      }
    end

    def restrict_constraints
      return {1 => 1} unless @restrict && model_restricted_attribute

      Sequel.negate(
        Sequel.qualify(@outgoing_alias, model_restricted_attribute.column_name) => true
      )
    end

    def model_restricted_attribute
      @outgoing_model.attributes[:restricted]
    end
  end
end
