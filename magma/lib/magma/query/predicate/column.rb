class Magma
  class ColumnPredicate < Magma::Predicate
    # This Predicate returns an actual attribute value of some kind - a number, integer, etc.,
    # or else a test on that value (number > 2, etc.)
    def initialize question, model, alias_name, attribute, parent_alias_name, is_subselect, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      @parent_alias_name = parent_alias_name
      @attribute = attribute
      @attribute_name = attribute.attribute_name.to_sym
      @column_name = attribute.column_name.to_sym
      @is_subselect = is_subselect
      process_args(query_args)
    end

    def self.inherited(subclass)
      Magma::Predicate.inherited(subclass)
    end

    def extract table, identity
      if @verb && @verb.gives?(:extract)
        @verb.do(:extract, table, identity)
      else
        table.first[column_name]
      end
    end

    def format
      default_format
    end

    def select
      @arguments.empty? ?
        # @is_subselect ?
        # [ Magma::TerminalSubselect.new(**subselect_params).coalesce ] :
        [ Sequel[alias_name][attribute_column_name].as(column_name) ] :
      []
    end

    def attribute_column_name
      @column_name
    end

    protected

    def column_name
      :"#{alias_name}_#{attribute_column_name}"
    end

    private

    def subselect_params
      require 'pry'
      binding.pry
      {
        incoming_alias: @parent_alias_name,
        attribute_alias: column_name,
        outgoing_alias: alias_name,
        outgoing_identifier_column_name: @model.identity.column_name,
        outgoing_model: @model,
        outgoing_column_name: attribute_column_name,
        restrict: @question.restrict?
      }
    end
  end
end
