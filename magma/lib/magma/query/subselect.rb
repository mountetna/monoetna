class Magma
  class SubselectError < Exception
  end

  class Subselect
    attr_reader :attribute_alias

    def initialize(
      parent_alias:,
      attribute_alias:,
      child_model:,
      child_alias:,
      child_identifier_column_name:,
      child_fk_column_name:,
      child_column_name: nil,
      subselect: nil
    )
      raise SubselectError.new("Cannot set both child_column_name and subselect.") unless child_column_name.nil? || subselect.nil?

      @parent_alias = parent_alias.to_sym
      @attribute_alias = attribute_alias.to_sym
      @child_model = child_model
      @child_alias = child_alias.to_sym
      @child_identifier_column_name = child_identifier_column_name.to_sym
      @child_fk_column_name = child_fk_column_name.to_sym
      @child_column_name = child_column_name
      @subselect = subselect
    end

    def coalesce
      require 'pry'
      binding.pry
      Sequel.function(
        :coalesce,
        inner_select,
        default_value
      ).as(attribute_alias)
    end

    def hash
      attribute_alias.hash
    end

    def eql?(other)
      attribute_alias == other.attribute_alias
    end

    private

    def inner_select
      subselect_query = @child_model.from(
        Sequel.as(@child_model.table_name, @child_alias)
      )

      subselect_query.select(
        Sequel.function(
          :json_agg,
          identifier_tuple
        )
      ).where(**subselect_constraint)
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
        Sequel.qualify(@child_alias, @child_fk_column_name) => Sequel.qualify(@parent_alias, :id)
      }
    end
  end
end
