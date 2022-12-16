class Magma
  class ColumnPredicate < Magma::Predicate
    # This Predicate returns an actual attribute value of some kind - a number, integer, etc.,
    # or else a test on that value (number > 2, etc.)
    def initialize question, model, alias_name, attribute, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      @attribute_obj = attribute
      @attribute_name = attribute.attribute_name.to_sym
      @column_name = attribute.column_name.to_sym
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
      require 'pry'
      binding.pry
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

    def root_model_column?
      @question.model == @model
    end

    def aggregate_function
      is_json_column? ? :string_agg : :json_agg
    end

    def is_collection_column?
      collection_column_types = [
        Magma::ChildAttribute,
        Magma::ForeignKeyAttribute
      ]

      collection_column_types.include?(@attribute_obj.class)
    end

    def is_json_column?
      json_column_types = [
        Magma::FileAttribute,
        Magma::ImageAttribute,
        Magma::FileCollectionAttribute
      ]

      json_column_types.include?(@attribute_obj.class)
    end

    def is_aggregated_attribute?
      !root_model_column? || is_json_column? || is_collection_column?
    end

    def aggregated_column_name
      "#{alias_name}_#{@attribute_name}_aggregated"
    end
  end
end
