class Magma
  class SubselectError < Exception
  end

  class Subselect
    def initialize(
      incoming_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      restrict:,
      outgoing_fk_column_name:,
      requested_data:,
      filters:,
      alias_column: true
    )
      @incoming_alias = incoming_alias.to_sym
      @outgoing_model = outgoing_model
      @outgoing_alias = outgoing_alias.to_sym
      @outgoing_identifier_column_name = outgoing_identifier_column_name.to_sym
      @restrict = restrict
      @outgoing_fk_column_name = outgoing_fk_column_name
      @requested_data = requested_data # another subselect, or a column mame on the outgoing model
      @filters = filters
      @alias_column = alias_column
    end

    def coalesce
      build_coalesce(inner_select)
    end

    def hash
      subselect_column_alias.hash
    end

    def eql?(other)
      subselect_column_alias == other.subselect_column_alias
    end

    def subselect_column_alias
      @requested_data.subselect_column_alias
    end

    private

    def build_coalesce(select_statement)
      main_data = Sequel.function(
        :coalesce,
        select_statement,
        default_value
      )

      main_data = main_data.as(subselect_column_alias) if @alias_column

      main_data
    end

    def inner_select
      query = subselect_query

      query = query.select(
        Sequel.function(
          :json_agg,
          subselect_data
        )
      ).where(
        **subselect_constraints
      ).where(
        restrict_constraints
      )

      apply_filter_constraints(query)
    end

    def subselect_query
      @outgoing_model.from(
        Sequel.as(@outgoing_model.table_name, @outgoing_alias)
      )
    end

    def subselect_data
      identifier_tuple
    end

    def identifier_tuple
      # Returns a tuple like [ ::identifier, <requested data> ]
      Sequel.function(:json_build_array,
        outgoing_identifier_column,
        outgoing_data
      )
    end

    def outgoing_identifier_column
      Sequel.qualify(@outgoing_alias, @outgoing_identifier_column_name)
    end

    def outgoing_data
      @requested_data.is_a?(Magma::Subselect) ?
        @requested_data.coalesce :
        Sequel.qualify(@outgoing_alias, @requested_data.outgoing_data_column_name)
    end

    def default_value
      Sequel.cast('[]', :json)
    end

    def subselect_constraints
      {
        Sequel.qualify(@outgoing_alias, @outgoing_fk_column_name) => Sequel.qualify(@incoming_alias, :id)
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

    def apply_filter_constraints(query)
      return query if @filters.empty?

      @filters.map do |filter|
        filter.flatten.map(&:constraint).inject(:+)
      end.inject(:+).each do |constraint|
        query = constraint.apply(query)
      end

      query
    end
  end
end
