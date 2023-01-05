class Magma
  class ColumnPredicate < Magma::Predicate
    # This Predicate returns an actual attribute value of some kind - a number, integer, etc.,
    # or else a test on that value (number > 2, etc.)
    def initialize question, model, alias_name, attribute, is_subselect, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
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
        @is_subselect ?
        [ Magma::Subselect.new(*subselect_params) ] :
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
      {
        parent_alias: 'how do you figure this out?',
        parent_id_column_name: parent_id_column,
        attribute_alias: column_name,
        child_alias: alias_name,
        child_identifier_column_name: @model.identity.column_name,
        child_fk_column_name: fk_attribute.column_name,
        child_table_name: @model.table_name,
        child_column_name: attribute_column_name
      }
    end

    def fk_attribute
      @fk_attribute ||= @model.attributes.values.select do |attribute|
        attribute.is_a?(Magma::ParentAttribute)
      end.first
    end

    def parent_id_column
      fk_attribute.link_model.attributes[fk_attribute.link_attribute_name.to_sym].column_name.to_sym
    end
  end
end
