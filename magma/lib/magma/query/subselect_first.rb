require_relative 'subselect'

class Magma
  class SubselectFirst < Magma::Subselect
    def coalesce
      limited_inner_select.as(subselect_column_alias)
    end

    private

    def limited_inner_select
      query = subselect_query

      query = query.select(
        outgoing_data
      ).where(
        **subselect_constraints
      ).where(
        restrict_constraints
      ).order_by(
        outgoing_table_id
      ).limit(1)

      apply_filter_constraints(query)
    end

    def outgoing_table_id
      Sequel.qualify(@outgoing_alias, :id)
    end
  end
end