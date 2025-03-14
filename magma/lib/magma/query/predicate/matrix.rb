require 'set'

class Magma
  class MatrixPredicate < Magma::Predicate
    attr_reader :attribute
    def initialize question, model, alias_name, attribute, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      @attribute = attribute
      @attribute_name = attribute.name.to_sym
      @column_name = attribute.column_name.to_sym
      @requested_identifiers = Set.new
      process_args(query_args)
    end

    verb '::slice', Array do
      child Array

      extract do |table, identity|
        answer = AnswerTuple.new(
          table.first[identity],
          table.first[column_name]
        )
        Magma::Answer.new(MatrixValue.new(self, answer.data, @arguments[1]))
      end
      validate do |_, validation_list|
        (validation_list - @predicate.attribute.validation_object.options).empty? && !validation_list.empty?
      end
      format { [ default_format, @arguments[1] ] }
    end

    verb nil do
      child Array

      extract do |table, identity|
        answer = AnswerTuple.new(
          table.first[identity],
          table.first[column_name]
        )
        Magma::Answer.new(MatrixValue.new(self, answer.data))
      end
      format { [ default_format, @attribute.validation_object.options ] }
    end

    def select(incoming_alias_name=nil, incoming_attribute=nil)
      [ select_column_name.as(column_name) ]
    end

    def matrix_row(data, column_names)
      @attribute.matrix_row_json(data, column_names)
    end

    def attribute_column_name
      @column_name
    end

    protected

    def column_name
      :"#{alias_name}_#{@column_name}"
    end

    def select_column_name
      Sequel[alias_name][@column_name.to_sym]
    end

    class MatrixValue
      def initialize(predicate, data, column_names=nil)
        @predicate = predicate
        @data = data
        @column_names = column_names
      end

      def to_json(options={})
        @predicate.matrix_row(@data,@column_names)
      end
    end
  end
end
