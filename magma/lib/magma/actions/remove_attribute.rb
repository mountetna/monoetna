class Magma
  class RemoveAttributeAction < BaseAction
    def perform
      model.attributes.delete(attribute.attribute_name.to_sym)
      attribute.delete
      @errors.empty?
    end

    def target_models
      if model
        [model]
      else
        []
      end
    end

    def validations
      [
        :validate_attribute_exists,
        :validate_not_primary_key,
        :validate_not_link,
        :validate_not_break_graph
      ]
    end

    def validate_not_primary_key
      return unless attribute
      return unless attribute.primary_key? || attribute.is_a?(Magma::IdentifierAttribute)

      @errors << Magma::ActionError.new(
        message: 'Cannot remove identity attribute',
        source: @action_params.slice(:attribute_name, :model_name)
      )
    end

    def validate_attribute_exists
      return if attribute

      @errors << Magma::ActionError.new(
        message: 'Attribute does not exist',
        source: @action_params.slice(:attribute_name, :model_name)
      )
    end

    def validate_not_link
      return unless attribute

      link_attributes = [
        Magma::LinkAttribute,
        Magma::ChildAttribute,
        Magma::CollectionAttribute
      ]

      return unless link_attributes.include?(attribute.class)

      # parent -> collection should raise a break-graph error, instead
      return if (attribute.is_a?(Magma::CollectionAttribute) || attribute.is_a?(Magma::ChildAttribute)) && attribute.link_model.attributes[attribute.link_attribute_name.to_sym].is_a?(Magma::ParentAttribute)

      @errors << Magma::ActionError.new(
        message: 'Use remove_link action to remove a link attribute',
        source: @action_params.slice(:attribute_name, :model_name)
      )
    end

    def validate_not_break_graph
      return unless attribute

      breaking_attributes = [
        Magma::ForeignKeyAttribute,
        Magma::TableAttribute,
        Magma::CollectionAttribute,
        Magma::ParentAttribute
      ]
      return unless breaking_attributes.include?(attribute.class)

      @errors << Magma::ActionError.new(
        message: 'Use reparent_model action to change this attribute',
        source: @action_params.slice(:attribute_name, :model_name)
      )
    end

    def model
      @model ||= Magma.instance.get_model(
        @project_name,
        @action_params[:model_name]
      )
    end

    def attribute
      return nil unless @action_params[:attribute_name].respond_to?(:to_sym)
      @attribute ||= model.attributes[@action_params[:attribute_name].to_sym]
    end
  end

  def model
    @model ||= Magma.instance.get_model(
      @project_name,
      @action_params[:model_name]
    )
  end

  def attribute
    return nil unless model
    return nil unless @action_params[:attribute_name].respond_to?(:to_sym)
    @attribute ||= model.attributes[@action_params[:attribute_name].to_sym]
  end
end
