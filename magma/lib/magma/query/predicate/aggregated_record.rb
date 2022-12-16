class Magma
  class AggregatedRecordPredicate < Magma::Predicate
    attr_reader :model

    def initialize question, model, alias_name, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      require 'pry'
      binding.pry
      process_args(query_args)
    end

    verb '::identifier' do
      child :record_child

      select_columns do
        [ column_name ]
      end
    end

    verb :attribute_name do
      child :record_child

      join :attribute_join

      select_columns do
        [ column_name ]
      end
    end

    def to_hash
      super.merge(
        model: model
      )
    end

    def select
      [
        Sequel.function(:json_agg, Sequel[@alias_name][column_name(@arguments[0])]).distinct.as(column_name_alias)
      ]
    end

    private

    def record_child
      RecordPredicate.new(@question, @model, @alias_name, *@query_args)
    end

    def column_name(attribute = @model.identity)
      if attribute.is_a?(String) || attribute.is_a?(Symbol)
        attribute = @model.attributes[attribute.to_sym]
        if attribute.nil?
          attribute = @model.identity
        end
      end

      attribute.column_name.to_sym
    end

    def column_name_alias
      "#{@alias_name}_#{column_name(@arguments[0])}_aggregated"
    end
  end
end
