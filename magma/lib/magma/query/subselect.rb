require_relative 'subselect_base'

class Magma
  class Subselect < Magma::SubselectBase
    def initialize(
      incoming_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      restrict:,
      outgoing_fk_column_name:,
      requested_data:,
      filters:,
      alias_column: true
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
      @alias_column = alias_column
    end

    def coalesce
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
      main_data = Sequel.function(
        :coalesce,
        select_statement,
        default_value
      )

      main_data = main_data.as(subselect_column_alias) if @alias_column

      main_data
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
      @requested_data.is_a?(Magma::Subselect) ?
        @requested_data.coalesce :
        Sequel.qualify(@outgoing_alias, @requested_data.outgoing_data_column_name)
    end

    def default_value
      Sequel.cast('[]', :json)
    end
  end
end
