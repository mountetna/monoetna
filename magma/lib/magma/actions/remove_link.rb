class Magma
  class RemoveLinkAction < BaseAction
    def perform
      reciprocal_model.attributes.delete(reciprocal_attribute.attribute_name.to_sym)
      model.attributes.delete(attribute.attribute_name.to_sym)
      attribute.delete
      reciprocal_attribute.delete
      @errors.empty?
    end

    def target_models
      if model
        [model, attribute.link_model]
      else
        []
      end
    end

    def validations
      [
        :validate_attribute_exists,
        :validate_is_link,
      ]
    end

    def validate_attribute_exists
      return if attribute

      @errors << Magma::ActionError.new(
        message: 'Attribute does not exist',
        source: @action_params.slice(:attribute_name, :model_name)
      )
    end

    def validate_is_link
      return unless attribute
      return if attribute.is_a?(Magma::LinkAttribute) || reciprocal_attribute&.is_a?(Magma::LinkAttribute)

      @errors << Magma::ActionError.new(
        message: 'Attribute is not part of a link',
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
      return nil unless model
      return nil unless @action_params[:attribute_name].respond_to?(:to_sym)
      @attribute ||= model.attributes[@action_params[:attribute_name].to_sym]
    end

    def reciprocal_model
      return nil unless attribute
      return nil unless attribute.link_model_name

      @reciprocal_model ||= Magma.instance.get_model(
        @project_name,
        attribute.link_model_name
      )
    end

    def reciprocal_attribute
      return nil unless reciprocal_model
      return nil unless attribute.link_attribute_name

      @reciprocal_attribute ||= reciprocal_model.attributes[attribute.link_attribute_name.to_sym]
    end
  end
end
