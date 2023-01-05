class Magma
  class SubselectError < Exception
  end

  class Subselect
    attr_reader :attribute_alias

    def initialize(
      parent_alias:,
      parent_id_column_name:,
      attribute_alias:,
      child_alias:,
      child_identifier_column_name:,
      child_fk_column_name:,
      child_table_name:,
      child_column_name: nil,
      subselect: nil
    )
      raise SubselectError.new("Cannot set both child_column_name and subselect.") unless child_column_name.nil? || subselect.nil?

      @parent_alias = parent_alias.to_sym
      @parent_id_column_name = parent_id_column_name.to_sym
      @attribute_alias = attribute_alias.to_sym
      @child_alias = child_alias.to_sym
      @child_identifier_column_name = child_identifier_column_name.to_sym
      @child_fk_column_name = child_fk_column_name.to_sym
      @child_table_name = child_table_name
      @child_column_name = child_column_name
      @subselect = subselect
    end

    def apply query
      require 'pry'
      binding.pry
      query.select_append(
        coalesce(
          inner_select,
          default_value
        ).as(attribute_alias)
      )
    end

    def hash
      attribute_alias.hash
    end

    def eql?(other)
      attribute_alias == other.attribute_alias
    end

    private

    def inner_select
      subselect_query = child_table.from(
        Sequel.as(@child_table_name, @child_alias)
      )

      subselect_query.select(
        Sequel.function(
          :json_agg,
          identifier_tuple
        )
      ).where(*subselect_constraint)

      subselect_query
    end

    def identifier_tuple
      # Returns a tuple like [ ::identifier, <requested value> ]
      Sequel.function(:json_build_array,
        child_identifier_column,
        @subselect || child_column
      )
    end

    def child_identifier_column
      Sequel.qualify(@child_alias, @child_identifier_column_name)
    end

    def child_column
      Sequel.qualify(@child_alias, @child_column_name)
    end

    def default_value
      Sequel.cast('[]', :json)
    end

    def subselect_constraint
      {
        Sequel.qualify(@child_alias, @child_fk_column_name) => Sequel.qualify(@parent_alias, @parent_id_column_name)
      }
    end
  end
end
