require_relative 'with_aggregated_data_module'

class Magma
  class ColumnPredicate < Magma::Predicate
    include WithAggregatedDataModule

    # This Predicate returns an actual attribute value of some kind - a number, integer, etc.,
    # or else a test on that value (number > 2, etc.)
    def initialize question, model, alias_name, attribute, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      @attribute = attribute
      @attribute_name = attribute.attribute_name.to_sym
      @column_name = attribute.column_name.to_sym
      process_args(query_args)
    end

    def self.inherited(subclass)
      Magma::Predicate.inherited(subclass)
    end

    def extract table, identity, is_all
      if @verb && @verb.gives?(:extract)
        @verb.do(:extract, table, identity)
      elsif table.first[column_name].nil?
        Magma::NilAnswer.new
      else
        Magma::Answer.new(table.first[column_name])
      end
    end

    def format
      default_format
    end

    def select
      @arguments.empty? ?
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
  end
end
