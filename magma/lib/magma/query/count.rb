require_relative 'with_restrict_module'
require_relative 'with_filters_module'
require_relative 'with_subqueries_module'

class Magma
  class Count
    include WithRestrictModule
    include WithFiltersModule
    include WithSubqueriesModule

    def initialize(model:, filters:, restrict:, table_alias_name:, subqueries:)
      @model = model
      @filters = filters
      @restict = restrict
      @table_alias_name = table_alias_name
      @subqueries = subqueries
    end

    def build
      count_select.as(count_column_alias)
    end

    def count_column_alias
      :"#{@table_alias_name}_count"
    end

    private

    def base_query
      @model.from(
        Sequel.as(@model.table_name, @table_alias_name)
      )
    end

    def count_select
      query = base_query

      query = query.select(
        count_column
      ).where(
        restrict_constraints
      )

      query = apply_filters(query)
      apply_subqueries(query)
    end

    def count_column
      Sequel.function(:count,
        Sequel.function(:distinct,
          Sequel.qualify(
            @table_alias_name,
            outgoing_column_name
          )
        )
      )
    end

    def restrict_constraints
      restrict_constraints_for_model_alias(@model, @table_alias_name)
    end

    def outgoing_column_name
      :id
    end
  end
end