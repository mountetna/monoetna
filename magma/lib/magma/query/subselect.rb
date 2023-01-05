class Magma
  class SubselectError < Exception
  end

  class Subselect
    attr_reader :attribute_alias

    def initialize(
      incoming_alias:,
      attribute_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      restrict:,
      outgoing_column_name:,
      subselect: nil
    )
      raise SubselectError.new("Must provide one of outgoing_column_name or subselect.") if outgoing_column_name.nil? && subselect.nil?
      require 'pry'
      binding.pry if incoming_alias.nil?
      @incoming_alias = incoming_alias.to_sym
      @attribute_alias = attribute_alias.to_sym
      @outgoing_model = outgoing_model
      @outgoing_alias = outgoing_alias.to_sym
      @outgoing_identifier_column_name = outgoing_identifier_column_name.to_sym
      @restrict = restrict
      @outgoing_column_name = outgoing_column_name
    end

    def coalesce
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
      subselect_query = @outgoing_model.from(
        Sequel.as(@outgoing_model.table_name, @outgoing_alias)
      )

      subselect_query.select(
        Sequel.function(
          :json_agg,
          subselect_data
        )
      ).where(**subselect_constraints).where(restrict_constraints)
    end

    def subselect_data
      raise SubselectError.new("Must be implemented in a subclass.")
    end

    def outgoing_identifier_column
      Sequel.qualify(@outgoing_alias, @outgoing_identifier_column_name)
    end

    def outgoing_column
      Sequel.qualify(@outgoing_alias, @outgoing_column_name)
    end

    def default_value
      Sequel.cast('[]', :json)
    end

    def subselect_constraints
      {
        Sequel.qualify(@outgoing_alias, @outgoing_column_name) => Sequel.qualify(@incoming_alias, :id)
      }
    end

    def restrict_constraints
      return {1 => 1} unless @restrict && model_restricted_attribute

      Sequel.negate(
        Sequel.qualify(@outgoing_alias, model_restricted_attribute.column_name) => true
      )
    end

    def model_restricted_attribute
      @outgoing_model.attributes[:restricted]
    end
  end
end
