require_relative 'subselect_base'
require_relative 'with_subqueries_module'

class Magma
  class SubselectCountBuilder < Magma::SubselectBase
    include WithSubqueriesModule

    def initialize(
      incoming_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      restrict:,
      outgoing_fk_column_name:,
      subqueries:,
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
      @subqueries = subqueries
    end

    def build
      count_select
    end

    def subselect_column_alias
      :"#{@outgoing_alias}_count"
    end

    private

    def count_select
      query = subselect_query

      query = query.select(
        count_column
      ).where(
        **subselect_constraints
      ).where(
        restrict_constraints
      )

    query = apply_filters(query)
    apply_subqueries(query)
    end

    def count_column
      Sequel.function(:count, outgoing_identifier_column)
    end
  end
end