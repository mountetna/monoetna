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
      [:validate_attribute_exists, :validate_not_primary_key]
    end

    def validate_not_primary_key
      return unless attribute
      p attribute
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

    def validate_restricted_attribute
      return if @action_params[:attribute_name] != 'restricted' || !@action_params[:restricted]

      @errors << Magma::ActionError.new(
        message: "restricted column may not, itself, be restricted",
        source: @action_params.slice(:project_name, :model_name, :attribute_name, :restricted)
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
