require_relative '../model'

class Magma
  class ModelSubselectPredicate < Magma::ModelPredicate
    def self.verbs
      Magma::ModelPredicate.verbs.merge(@verbs)
    end

    verb '::first' do
      child :record_child
      extract do |table,return_identity|
        if table.empty?
          Magma::NilAnswer.new
        else # there is only one row in the table
          child_extract(table, identity)
        end
      end
      format do
        child_format
      end

      select_columns do |incoming_alias_name, incoming_attribute|
        [
          Magma::SubselectFirstBuilder.new(**subselect_params(
            incoming_alias_name,
            incoming_attribute
          ))
        ]
      end
    end

    verb '::all' do
      child :record_child
      extract do |table,return_identity|
        root_table_identifier = table.first[return_identity]
        child_answer = child_predicate.extract(table, identity, true)

        if root_table_identifier.nil?
          child_answer
        else
          Magma::AnswerTuple.new(
            root_table_identifier,
            child_answer
          )
        end
      end
      format do
        [
          default_format,
          child_format
        ]
      end
    end

    verb '::count' do
      child Numeric
      extract do |table,return_identity|
        Magma::Answer.new(table.first[count_column_name])
      end
      format { 'Numeric' }

      select_columns do |incoming_alias_name=nil, incoming_attribute=nil|
        if incoming_alias_name && incoming_attribute
          [
            Magma::SubselectCountBuilder.new(**base_subselect_params(
              incoming_alias_name,
              incoming_attribute
            ))
          ]
        else
          []
        end
      end
    end

    def record_child
      RecordSubselectPredicate.new(@question, @model, alias_name, *@query_args)
    end

    def select(incoming_alias_name, incoming_attribute)
      if @verb && @verb.gives?(:select_columns)
        @verb.do(:select_columns, incoming_alias_name, incoming_attribute)
      else
        [
          Magma::SubselectBuilder.new(**subselect_params(
            incoming_alias_name,
            incoming_attribute
          ))
        ]
      end
    end

    private

    def child_subselect(incoming_attribute)
      child_predicate.select(alias_name, incoming_attribute).first
    end

    def subselect_params(incoming_alias_name, incoming_attribute)
      base_subselect_params(incoming_alias_name, incoming_attribute).update({
        requested_data: child_subselect(child_predicate_incoming_attribute(incoming_attribute))
      })
    end

    def base_subselect_params(incoming_alias_name, incoming_attribute)
      {
        incoming_alias: incoming_alias_name,
        outgoing_alias: alias_name,
        outgoing_identifier_column_name: @model.identity.column_name,
        outgoing_fk_column_name: outgoing_attribute(incoming_attribute).column_name,
        outgoing_model: @model,
        restrict: @question.restrict?,
        filters: @filters,
      }
    end

    def child_predicate_incoming_attribute(incoming_attribute)
      child_predicate.attribute || incoming_attribute
    end

    def outgoing_attribute(incoming_attribute)
      incoming_attribute.link_model.attributes[incoming_attribute.link_attribute_name.to_sym]
    end
  end
end