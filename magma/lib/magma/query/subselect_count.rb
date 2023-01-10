require_relative 'subselect_base'

class Magma
  class SubselectCount < Magma::SubselectBase

    def build
      count_select.as(count_column_alias)
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

      apply_filter_constraints(query)
    end

    def count_column
      Sequel.function(:count, outgoing_identifier_column)
    end

    def count_column_alias
      :"#{@outgoing_alias}_count"
    end
  end
end