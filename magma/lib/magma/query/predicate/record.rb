class Magma
  class RecordPredicate < Magma::Predicate
  # This object takes several arguments:
  #   1) It can accept any of its attributes as arguments
  #      Here are the Magma attribute types:
  #        ChildAttribute - this returns another Record predicate
  #        CollectionAttribute
  #        TableAttribute - these both return a Model predicate
  #        FileAttribute - this returns a File predicate
  #        ImageAttribute - this returns a Image predicate
  #        ForeignKey - this returns a Record predicate
  #        Attribute - this, depending on its type, can have different results
  #          If the type is a String, you get a String predicate
  #          If the type is an Integer or Float you get a Number predicate
  #          if the type is a DateTime you get a DateTime predicate
  #          if the type is a Boolean you get a Boolean predicate
  #   2) ::has
  #   3) ::identifier
    attr_reader :model

    def initialize question, model, alias_name, *query_args
      super(question)
      @model = model
      @alias_name = alias_name
      process_args(query_args)
    end

    verb '::identifier' do
      child do
        attribute_child(@model.identity.attribute_name)
      end

      join :attribute_join

      select_columns :attribute_select
    end

    verb '::has', :attribute_name do
      child TrueClass

      constraint do
        attribute = valid_attribute(@arguments[1])
        case attribute
        when Magma::ForeignKeyAttribute, Magma::LinkAttribute
          not_null_constraint(attribute.foreign_id)
        when Magma::TableAttribute, Magma::CollectionAttribute, Magma::ChildAttribute
          # do nothing in this case, instead, we'll add a join on the
          #   child table
        else
          not_null_constraint(attribute.column_name.to_sym)
        end
      end

      join do
        attribute = valid_attribute(@arguments[1])
        case attribute
        when Magma::TableAttribute, Magma::CollectionAttribute, Magma::ChildAttribute
          inner_join_child(attribute)
        end
      end
    end

    verb '::lacks', :attribute_name do
      child TrueClass

      constraint do
        attribute = valid_attribute(@arguments[1])
        case attribute
        when Magma::ForeignKeyAttribute, Magma::LinkAttribute
          null_constraint(attribute.foreign_id)
        when Magma::FileAttribute, Magma::ImageAttribute
          or_constraint([
            json_constraint(attribute.column_name.to_sym, "filename", nil),
            # The string "null" appears if someone fetches a ::temp URL
            #   for a file attribute, because the FileSerializer returns
            #   nil to the loader. This gets translated into JSON "null"
            #   and saved to the database, instead of a SQL NULL.
            # So when we check for ::lacks here, we also need to check for
            #   the string "null".
            json_constraint(attribute.column_name.to_sym, "filename", "null"),
            null_constraint(attribute.column_name.to_sym),
          ])
        when Magma::TableAttribute, Magma::CollectionAttribute, Magma::ChildAttribute
          not_constraint(
            :id,
            lacks_child_subquery(attribute)
          )
        else
          null_constraint(attribute.column_name.to_sym)
        end
      end
    end

    verb '::metrics' do
      child do
        Magma::MetricsPredicate.new(@question, @model, alias_name, *@query_args)
      end
    end

    verb :attribute_name do
      child do
        attribute_child(@arguments[0])
      end
      join :attribute_join

      select_columns :attribute_select
    end

    verb Array do
      child do
        Magma::TablePredicate.new(@question, @model, alias_name, @arguments[0], *@query_args)
      end
    end

    def to_hash
      super.merge(
        model: model
      )
    end

    def select(incoming_alias_name=nil, incoming_attribute=nil)
      if @verb && @verb.gives?(:select_columns)
        @verb.do(:select_columns, incoming_alias_name, incoming_attribute)
      else
        super()
      end
    end

    def attribute
      return nil if @arguments.empty?

      valid_attribute(@arguments[0])
    end

    def collection_attribute?
      child_predicate.is_a?(Magma::FileCollectionPredicate)
    end

    private

    def attribute_select(incoming_alias_name=nil, incoming_attribute=nil)
      [
        child_predicate.select(
          alias_name,
          valid_attribute(@arguments[0])
        )
      ].flatten
    end

    def attribute_name(argument)
      @model.has_attribute?(argument) || argument == :id || argument == '::identifier'
    end

    def attribute_alias
      # We don't have a predicates for collection-type attributes,
      #   so we provide a convenience method to generate new aliases
      10.times.map{ (97+rand(26)).chr }.join.to_sym
    end

    def inner_join_child(attribute)
      Magma::Join.new(
        # left table
        table_name,
        alias_name,
        :id,

        #right table
        attribute.link_model.table_name,
        attribute_alias,
        collection_attribute_fk_column(attribute),
        inner_join: true
      )
    end

    def collection_attribute_fk_column(attribute)
      attribute.link_model.attributes[attribute.link_attribute_name.to_sym].column_name.to_sym
    end

    def lacks_child_subquery(attribute)
      attribute_model_alias = attribute_alias
      fk_column = collection_attribute_fk_column(attribute)
      attribute.link_model.select(
        Sequel.as(
          Sequel.qualify(attribute_model_alias, fk_column),
          fk_column
        )
      ).from(
        Sequel.as(attribute.link_model.table_name, attribute_model_alias)
      ).where(
        Sequel.negate(
          Sequel.qualify(attribute_model_alias, fk_column) => nil
        )
      )
      # This "Where" clause is important to apply because any NULL fk records
      #   on the child table will cause the subquery result to include NULL,
      #   which means the "NOT IN" constraint ambiguity causes the overall
      #   query to return []. https://stackoverflow.com/q/1406215
    end

    def attribute_join
      attribute = valid_attribute(@arguments[0])
      case attribute
      when Magma::ForeignKeyAttribute
        return Magma::Join.new(
          # left table
          table_name,
          alias_name,
          attribute.foreign_id,

          #right table
          @child_predicate.table_name,
          @child_predicate.alias_name,
          :id
        )
      when Magma::ChildAttribute
        return Magma::Join.new(
          #left table
          table_name,
          alias_name,
          :id,

          #right table
          @child_predicate.table_name,
          @child_predicate.alias_name,
          attribute.foreign_id,
        )
      end
    end

    def attribute_child(attribute_name)
      attribute = valid_attribute(attribute_name)
      if @question.restrict? && attribute.respond_to?(:restricted) && attribute.restricted
        raise Etna::Forbidden, "Cannot query for restricted attribute #{attribute_name}"
      end
      case attribute
      when :id
        return Magma::NumberPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::ForeignKeyAttribute
        return Magma::RecordPredicate.new(@question, attribute.link_model, nil, *@query_args)
      when Magma::ChildAttribute
        return Magma::ChildModelPredicate.new(@question, attribute.link_model, nil, *@query_args)
      when Magma::TableAttribute, Magma::CollectionAttribute
        return Magma::ModelSubselectPredicate.new(@question, attribute.link_model, *@query_args)
      when Magma::FileAttribute, Magma::ImageAttribute
        return Magma::FilePredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::FileCollectionAttribute
        return Magma::FileCollectionPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::MatchAttribute
        return Magma::MatchPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::MatrixAttribute
        return Magma::MatrixPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::StringAttribute
        return Magma::StringPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::IntegerAttribute, Magma::FloatAttribute
        return Magma::NumberPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::DateTimeAttribute, Magma::ShiftedDateTimeAttribute
        return Magma::DateTimePredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::BooleanAttribute
        return Magma::BooleanPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      else
        invalid_argument! attribute.name
      end
    end

    def valid_attribute(attribute_name)
      case attribute_name
      when :id
        return :id
      when '::identifier'
        return @model.identity
      else
        return @model.attributes[attribute_name.to_sym]
      end
    end
  end
end
