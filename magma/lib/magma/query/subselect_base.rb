require_relative 'with_restrict_module'
require_relative 'with_filter_constraints_module'

class Magma
  class SubselectError < Exception
  end

  class SubselectBase
    include WithRestrictModule
    include WithFilterConstraintsModule

    def initialize(
      incoming_alias:,
      outgoing_model:,
      outgoing_alias:,
      outgoing_identifier_column_name:,
      restrict:,
      outgoing_fk_column_name:,
      filters:
    )
      @incoming_alias = incoming_alias.to_sym
      @outgoing_model = outgoing_model
      @outgoing_alias = outgoing_alias.to_sym
      @outgoing_identifier_column_name = outgoing_identifier_column_name.to_sym
      @restrict = restrict
      @outgoing_fk_column_name = outgoing_fk_column_name
      @filters = filters
    end

    def hash
      subselect_column_alias.hash
    end

    def eql?(other)
      subselect_column_alias == other.subselect_column_alias
    end

    private

    def subselect_query
      @outgoing_model.from(
        Sequel.as(@outgoing_model.table_name, @outgoing_alias)
      )
    end

    def subselect_constraints
      {
        Sequel.qualify(@outgoing_alias, @outgoing_fk_column_name) => Sequel.qualify(@incoming_alias, :id)
      }
    end

    def restrict_constraints
      restrict_constraints_for_model_alias(@outgoing_model, @outgoing_alias)
    end

    def outgoing_identifier_column
      Sequel.qualify(@outgoing_alias, @outgoing_identifier_column_name)
    end

  end
end
