class Magma
  class ModelPredicate < Magma::Predicate
    # Model predicate - this is what the query will start with, probably
    #
    # "sample"
    #
    # This is a request for all objects of type "sample", so it's return type should be:
    #   [ Sample ]
    #

    # This object takes several arguments:
    #   1) It can accept an arbitrary list of filters, which are
    #      in the form of lists, e.g.:
    #
    #      [ "patient", "experiment", "name", "::equals", "Colorectal" ]
    #      [ "patient", "clinical", "parameter", [ "name", "::equals", "Gender" ], "::first", "value", "::equals", "Male" ]
    #
    #      Each one of these filters must reduce to a Boolean, or else it is
    #      invalid.  They must come first.
    #
    #   2) It can be reduced by a list operator. The list operators are:
    #      ::any - a Boolean that returns true if the list is non-zero
    #      ::first - returns the first item in the list, namely a Model
    #      ::all - returns every item in the list, represented by a Model
    #      ::count - returns the number of items in the list
    #      ::every - a Boolean that returns true if every item in the list is non-zero
    #      ::distinct - returns distinct values of items in the list. Only works with model -> attribute

    attr_reader :model

    def each_ancestor
      current_model = @model
      ancestral_path = []
      while current_model do
        yield current_model, ancestral_path
        ancestral_path << current_model.parent_model_name&.to_s
        current_model = current_model.parent_model
      end
    end

    def initialize(question, model, is_subselect, *query_args)
      super(question)
      @model = model
      @filters = []
      @is_subselect = is_subselect

      subquery_args, filter_args = Magma::SubqueryUtils.partition_args(self, query_args)

      subquery_args.each do |args|
        create_subquery(args)
      end

      # Since we are shifting off the the first elements on the query_args array
      # we look to see if the first element is an array itself. If it is then we
      # add it to the filters.
      while filter_args.first.is_a?(Array)
        create_filter(filter_args.shift)
      end

      add_filters

      process_args(filter_args)
    end

    verb '::first' do
      child :record_child
      extract do |table,return_identity|
        if table.empty?
          Magma::NilAnswer.new
        elsif @is_subselect # there is only one row in the table
          child_extract(table, identity)
        else
          child_extract(
            table.group_by do |row|
              row[identity]
            end.first&.last,
            identity
          )
        end
      end
      format do
        child_format
      end

      select_columns :select_first_column
    end

    verb '::all' do
      child :record_child
      extract do |table,return_identity|
        if @is_subselect # there is only one row in the table now
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
        else
          as_answer_tuple_array(
            table.group_by do |row|
              row[identity]
            end.map do |identifier,rows|
              next unless identifier
              [ identifier, child_extract(rows, identity, true) ]
            end.compact
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

    verb '::any' do
      child TrueClass

      subquery :join_subqueries

      subquery_config Magma::SubqueryConfig.new(magma_class: Magma::SubqueryInner, condition: "> 0")

      extract do |table,return_identity|
        Magma::Answer.new(
          table.any? do |row|
            row[identity]
          end
        )
      end
      format { 'Boolean' }
    end

    verb '::every' do
      child TrueClass

      subquery :join_subqueries

      subquery_config Magma::SubqueryConfig.new(magma_class: Magma::SubqueryInner, condition: "= count(*)")

      extract do |table,return_identity|
        Magma::Answer.new(
          table.length > 0 && table.all? do |row|
            row[identity]
          end
        )
      end
      format { 'Boolean' }
    end

    verb '::count' do
      child Numeric
      extract do |table,return_identity|
        Magma::Answer.new(table.first[count_column_name])
      end
      format { 'Numeric' }

      select_columns :select_count_column
    end

    verb '::distinct' do
      child :record_child

      constraint do
        validate_distinct

        distinct_constraint(record_child.child_predicate.attribute_column_name)
      end

      select_columns do
        []
      end

      extract do |table|
        Magma::Answer.new(
          table.map do |row| row.values end.flatten.uniq.compact
        )
      end
      format { [ child_format ] }
    end

    def create_filter(args)
      filter = FilterPredicate.new(
        question: @question,
        model: @model,
        alias_name: alias_name,
        query_args: args)

      unless filter.reduced_type == TrueClass
        raise ArgumentError,
          "Filter #{filter} does not reduce to Boolean #{filter.argument} #{filter.reduced_type}!"
      end

      @filters.push(filter)
    end

    def create_subquery(args)
      verb = self.class.match_verbs(args, self, true).first
      subquery = SubqueryOperator.new(
        predicate: self,
        question: @question,
        model_alias_name: alias_name,
        query_args: args,
        subquery_class: verb.do(:subquery_config).magma_class
      )

      @subqueries.push(subquery)
    end

    def add_filters
      if @question.restrict? && !@is_subselect
        # the model can be restricted, and we should withhold restricted data
        # Subselects handle their own restrictions without filters
        each_ancestor do |restriction_model, ancestors|
          if restriction_model.has_attribute?(:restricted)
            create_filter(ancestors + [ 'restricted', '::untrue' ])
          end
        end
      end
    end

    def record_child
      RecordPredicate.new(@question, @model, alias_name, @is_subselect, *@query_args)
    end

    def join
      join_filters
    end

    def subquery
      join_subqueries.concat(join_filter_subqueries)
    end

    def select
      if @verb && @verb.gives?(:select_columns)
        @verb.do(:select_columns)
      else
        @is_subselect ? [] : [ column_name.as(identity) ]
      end
    end

    def column_name(attribute = @model.identity)
      if attribute.is_a?(String) || attribute.is_a?(Symbol)
        attribute = @model.attributes[attribute.to_sym]
        if attribute.nil?
          attribute = @model.identity
        end
      end

      Sequel[alias_name][attribute.column_name.to_sym]
    end

    def constraint
      filter_constraints.concat(verb_constraints)
    end

    def filter_constraints
      @filters.map do |filter|
        filter.flatten.map(&:constraint).inject(&:+) || []
      end.inject(&:+) || []
    end

    def verb_constraints
      if @verb && @verb.gives?(:constraint)
        [ @verb.do(:constraint) ].compact
      else
        []
      end
    end

    def to_hash
      super.merge(
        model: model,
        filters: @filters.map do |filter|
          filter.flatten.map do |pred|
            pred.to_hash
          end
        end
      )
    end

    def identity
      alias_for_column(@model.identity.column_name)
    end

    def alias_for_column(column_name)
      :"#{alias_name}_#{column_name}"
    end

    def alias_for_attribute(attr)
      if attr.is_a?(String) || attr.is_a?(Symbol)
        attr = @model.attributes[attr.to_sym]
        if attr.nil?
          return identity
        end
      end

      alias_for_column(attr.column_name)
    end

    def generate_subselect(incoming_alias_name, incoming_attribute)
      return child_predicate.select unless @is_subselect

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

    def as_answer_tuple_array(raw_answer_tuples)
      Magma::AnswerTupleArray.from_raw_answer_tuples(raw_answer_tuples)
    end

    def count_column_name
      :"#{alias_name}_count"
    end

    def select_first_column(incoming_alias_name=nil, incoming_attribute=nil)
      if @is_subselect && incoming_alias_name && incoming_attribute
        [
          Magma::SubselectFirstBuilder.new(**subselect_params(
            incoming_alias_name,
            incoming_attribute
          ))
        ]
      else
        [ ]
      end
    end

    def select_count_column(incoming_alias_name=nil, incoming_attribute=nil)
      if @is_subselect && incoming_alias_name && incoming_attribute
        [
          Magma::SubselectCountBuilder.new(**base_subselect_params(
            incoming_alias_name,
            incoming_attribute
          ))
        ]
      elsif @is_subselect
        []
      else
        [ Magma::Count.new(
          model: @model,
          filters: @filters,
          table_alias_name: alias_name,
          restrict: @question.restrict?
        ).build ]
      end
    end

    def child_subselect(incoming_attribute)
      child_predicate.generate_subselect(alias_name, incoming_attribute).first
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

    def validate_distinct
      raise ArgumentError.new("Can only use ::distinct on string attributes.") unless record_child.child_predicate.is_a?(Magma::StringPredicate)
    end
  end
end
