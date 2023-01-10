require_relative 'with_restrict_module'
require_relative 'with_filter_constraints_module'

class Magma
  class Count
    include WithRestrictModule
    include WithFilterConstraintsModule

    def initialize(model:, filters:, restrict:, table_alias_name:)
      @model = model
      @filters = filters
      @restict = restrict
      @table_alias_name = table_alias_name
    end

    def build
      count_select.as(count_column_alias)
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
      ).distinct

      apply_filter_constraints(query)
    end

    def count_column
      Sequel.function(:count, @model.identity.column_name)
    end

    def restrict_constraints
      restrict_constraints_for_model_alias(@model, @table_alias_name)
    end

    def count_column_alias
      :"#{@table_alias_name}_count"
    end
  end
end