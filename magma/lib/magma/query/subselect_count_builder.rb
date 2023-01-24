require_relative 'subselect_base'

class Magma
  class SubselectCountBuilder < Magma::SubselectBase

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

      apply_filters(query)
    end

    def count_column
      Sequel.function(:count, outgoing_identifier_column)
    end
  end
end