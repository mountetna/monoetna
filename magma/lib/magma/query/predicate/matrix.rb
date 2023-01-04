require 'set'

class Magma
  class MatrixPredicate < Magma::Predicate
    attr_reader :attribute
    def initialize question, model, alias_name, attribute, should_aggregate, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      @attribute = attribute
      @attribute_name = attribute.name.to_sym
      @column_name = attribute.column_name.to_sym
      # @requested_identifiers = Set.new
      @should_aggregate = should_aggregate
      process_args(query_args)
    end

    verb '::slice', Array do
      child Array

      extract do |table, identity|
        MatrixValue.new(self, table.first[column_name]&.first, @arguments[1])
      end
      validate do |_, validation_list|
        (validation_list - @predicate.attribute.validation_object.options).empty? && !validation_list.empty?
      end
      format { [ default_format, @arguments[1] ] }
    end

    verb nil do
      child Array

      extract do |table, identity|
        MatrixValue.new(self, table.first[column_name]&.first)
      end
      format { [ default_format, @attribute.validation_object.options ] }
    end

    def select
      @should_aggregate ?
        [ Sequel.function(:json_agg, Sequel[alias_name][@column_name]).distinct.as(column_name) ] :
        [ Sequel[alias_name][@column_name].as(column_name) ]
    end

    def matrix_row(identifier, column_names)
      @attribute.matrix_row_json(identifier, column_names)
    end

    protected

    def column_name
      :"#{alias_name}_#{@column_name}"
    end

    class MatrixValue
      def initialize(predicate, data, column_names=nil)
        @predicate = predicate
        @data = data
        @column_names = column_names
      end

      def to_json(options={})
        @predicate.matrix_row(@data, @column_names)
      end
    end
  end
end
