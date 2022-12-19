class Magma
  class ColumnPredicate < Magma::Predicate
    # This Predicate returns an actual attribute value of some kind - a number, integer, etc.,
    # or else a test on that value (number > 2, etc.)
    def initialize question, model, alias_name, attribute, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      @upstream_aggregation = false
      set_attribute(attribute)
      process_args(query_args)
    end

    def self.inherited(subclass)
      Magma::Predicate.inherited(subclass)
    end

    def self.from_record_predicate_and_attribute(record_predicate, attribute)
      case attribute
      when :id
        predicate_class = Magma::NumberPredicate
      when Magma::FileAttribute, Magma::ImageAttribute
        predicate_class = Magma::FilePredicate
      when Magma::FileCollectionAttribute
        predicate_class = Magma::FileCollectionPredicate
      when Magma::MatchAttribute
        predicate_class = Magma::MatchPredicate
      when Magma::MatrixAttribute
        predicate_class = Magma::MatrixPredicate
      when Magma::StringAttribute
        predicate_class = Magma::StringPredicate
      when Magma::IntegerAttribute, Magma::FloatAttribute
        predicate_class = Magma::NumberPredicate
      when Magma::DateTimeAttribute, Magma::ShiftedDateTimeAttribute
        predicate_class = Magma::DateTimePredicate
      when Magma::BooleanAttribute
        predicate_class = Magma::BooleanPredicate
      else
        invalid_argument! attribute.name
      end

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

    def set_upstream_aggregation(upstream_aggregation)
      @upstream_aggregation = upstream_aggregation
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
      binding.pry if @upstream_aggregation
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
        Magma::ForeignKeyAttribute,
        Magma::LinkAttribute
      ]

      collection_column_types.include?(@attribute_obj.class)
    end

    def is_json_column?
      # Postgres cannot group_by JSON columns, so we apply a
      #   string_agg on them, instead.
      json_column_types = [
        Magma::FileAttribute,
        Magma::ImageAttribute,
        Magma::FileCollectionAttribute
      ]

      json_column_types.include?(@attribute_obj.class)
    end

    def is_aggregated_attribute?
      !root_model_column? || is_json_column? || is_collection_column? || @upstream_aggregation
    end

    def aggregated_column_name
      "#{alias_name}_#{@attribute_name}_aggregated"
    end
  end
end
