class Magma
  class RecordSubselectPredicate < Magma::RecordPredicate
    def self.verbs
      Magma::RecordPredicate.verbs
    end

    def collection_attribute?
      child_predicate.is_a?(Magma::FileCollectionSubselectPredicate)
    end

    private

    def attribute_select(incoming_alias_name, incoming_attribute=nil)
      incoming_attribute = valid_attribute(@arguments[0]) if incoming_attribute.nil?
      [ child_predicate.select(incoming_alias_name, incoming_attribute) ].flatten
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
      end
    end

    def attribute_child(attribute_name)
      attribute = valid_attribute(attribute_name)
      if @question.restrict? && attribute.respond_to?(:restricted) && attribute.restricted
        raise Etna::Forbidden, "Cannot query for restricted attribute #{attribute_name}"
      end
      case attribute
      when :id
        return Magma::NumberSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::ForeignKeyAttribute
        return Magma::RecordSubselectPredicate.new(@question, attribute.link_model, nil, *@query_args)
      when Magma::ChildAttribute
        return Magma::ChildModelSubselectPredicate.new(@question, attribute.link_model, nil, *@query_args)
      when Magma::TableAttribute, Magma::CollectionAttribute
        return Magma::ModelSubselectPredicate.new(@question, attribute.link_model, *@query_args)
      when Magma::FileAttribute, Magma::ImageAttribute
        return Magma::FileSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::FileCollectionAttribute
        return Magma::FileCollectionSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::MatchAttribute
        return Magma::MatchSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::MatrixAttribute
        return Magma::MatrixSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::StringAttribute
        return Magma::StringSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::IntegerAttribute, Magma::FloatAttribute
        return Magma::NumberSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::DateTimeAttribute, Magma::ShiftedDateTimeAttribute
        return Magma::DateTimeSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      when Magma::BooleanAttribute
        return Magma::BooleanSubselectPredicate.new(@question, @model, alias_name, attribute, *@query_args)
      else
        invalid_argument! attribute.name
      end
    end
  end
end
