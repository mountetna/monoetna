class Magma
  class Subselect
    def initialize(parent_table, parent_table_alias, parent_table_id, child_table, child_table_alias, child_table_id)
      @child_table = child_table
      @child_table_alias = child_table_alias.to_sym
      @child_table_id = child_table_id.to_sym

      @parent_table = parent_table
      @parent_table_alias = parent_table_alias.to_sym
      @parent_table_id = parent_table_id.to_sym

      @constraints = [
        { parent_table_column => child_table_column }
      ]
    end

    def apply query
      require 'pry'
      binding.pry
      query.select_append(
        coalesce(
          Sequel.function(:json_agg,
            Sequel.function(:json_build_array,
              @child_table_id,
              '' # column_to_aggregate -- or another subselect?
          )).distinct,
          Sequel.cast('[]', :json)
        ).as(child_table_column)
      )
    end

    def child_table_column
      Sequel.qualify(@child_table_alias, @child_table_id)
    end

    def parent_table_column
      Sequel.qualify(@parent_table_alias, @parent_table_id)
    end

    def hash
      child_table_column.hash + parent_table_column.hash
    end

    def eql?(other)
      child_table_column == other.child_table_column && parent_table_column == other.parent_table_column
    end
  end
end
