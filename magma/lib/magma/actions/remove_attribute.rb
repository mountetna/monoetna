class Magma
  class UpdateAttributeAction < BaseAction
    def perform
      # attribute.update(@action_params.slice(*Magma::Attribute::EDITABLE_OPTIONS))
    rescue Sequel::ValidationFailed => e
      Magma.instance.logger.log_error(e)

      @errors << Magma::ActionError.new(
        message: 'Update attribute failed',
        source: @action_params.slice(:attribute_name, :model_name),
        reason: e
      )

      attribute.initial_values.keys.each { |name| attribute.reset_column(name) }
    ensure
      return @errors.empty?
    end

    def target_models
      if model
        [model]
      else
        []
      end
    end

    private

    def validations
      [
        :validate_params,
        :validate_attribute_exists,
        :validate_not_identifier,
        :validate_not_link,
      ]
    end

    def validate_params
      return if @action_params[:attribute_name] && @action_params[:model_name]

      @errors << Magma::ActionError.new(
        message: 'model_name and attribute_name are required, but missing.',
        source: @action_params.slice(:attribute_name, :model_name)
      )
    end

    def validate_not_identifier
      return unless attribute
      return unless attribute.is_a?(Magma::IdentifierAttribute)

      @errors << Magma::ActionError.new(
        message: 'Cannot remove an identifier attribute',
        source: @action_params.slice(:attribute_name, :model_name)
      )
    end

    def validate_not_link
      return unless attribute
      return unless attribute.is_a?(Magma::IdentifierAttribute)

      @errors << Magma::ActionError.new(
        message: 'Cannot remove an identifier attribute',
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

    def model
      @model ||= Magma.instance.get_model(
        @project_name,
        @action_params[:model_name]
      )
    end

    def attribute
      return nil if model.nil?
      return nil if @action_params[:attribute_name].nil?
      @attribute ||= model.attributes[@action_params[:attribute_name].to_sym]
    end
  end
end
