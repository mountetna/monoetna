class Magma
  class ColumnPredicate < Magma::Predicate
    # This Predicate returns an actual attribute value of some kind - a number, integer, etc.,
    # or else a test on that value (number > 2, etc.)
    attr_reader :upstream_aggregation

    def initialize question, model, alias_name, attribute, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      set_attribute(attribute)
      process_args(query_args)
    end

    def self.inherited(subclass)
      Magma::Predicate.inherited(subclass)
    end

    def self.from_record_predicate_and_attribute(record_predicate, attribute)
      predicate_class = :id == attribute ?
        Magma::NumberPredicate :
        attribute.predicate_class

      return predicate_class.new(
        record_predicate.question,
        record_predicate.model,
        record_predicate.alias_name,
        attribute,
        *record_predicate.query_args
      )
    end

    def set_attribute(attribute)
      @attribute_obj = attribute
      @attribute_name = attribute.attribute_name.to_sym
      @column_name = attribute.column_name.to_sym
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
      binding.pry if @attribute_name == "::identifier"
      @arguments.empty? ?
        is_aggregated_attribute? ?
          [ Sequel.function(aggregate_function, column_to_aggregate).distinct.as(aggregated_column_name) ] :
          [ Sequel[alias_name][@column_name].as(column_name) ] :

        []
    end

    def attribute_column_name
      @column_name
    end

    protected

    def column_name
      :"#{alias_name}_#{@column_name}"
    end

    private

    def column_to_aggregate
      is_json_column? ?
        Sequel[alias_name][@column_name].cast_string(:text) :
        Sequel[alias_name][@column_name]
    end
  end
end
