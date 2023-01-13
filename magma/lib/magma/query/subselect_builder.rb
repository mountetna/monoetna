require_relative 'subselect_base'

class Magma
  class SubselectBuilder < Magma::SubselectBase
    def initialize(
      incoming_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      restrict:,
      outgoing_fk_column_name:,
      requested_data:,
      filters:
    )
      super(
        incoming_alias: incoming_alias,
        outgoing_model: outgoing_model,
        outgoing_alias: outgoing_alias,
        outgoing_identifier_column_name: outgoing_identifier_column_name,
        restrict: restrict,
        outgoing_fk_column_name: outgoing_fk_column_name,
        filters: filters,
      )
      @requested_data = requested_data # another subselect, or a column mame on the outgoing model
    end

    def build
      build_coalesce(inner_select)
    end

    def subselect_column_alias
      @requested_data.subselect_column_alias
    end

    def outgoing_data_column_name
      @requested_data.outgoing_data_column_name
    end

    private

    def build_coalesce(select_statement)
      Sequel.function(
        :coalesce,
        select_statement,
        default_value
      )
    end

    def inner_select
      query = subselect_query

      query = query.select(
        inner_select_contents
      ).where(
        **subselect_constraints
      ).where(
        restrict_constraints
      )

      apply_filters(query)
    end

    def inner_select_contents
      Sequel.function(
        :json_agg,
        subselect_data
      )
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

    def outgoing_data
      @requested_data.is_a?(Magma::SubselectBuilder) ?
        @requested_data.build :
        Sequel.qualify(@outgoing_alias, @requested_data.outgoing_data_column_name)
    end

    def default_value
      Sequel.cast('[]', :json)
    end
  end
end
