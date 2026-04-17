require_relative "./with_update_model_module"
require_relative "./with_template_validation_module"

class Magma
  class SetModelTemplateAction < BaseAction
    include WithUpdateModel
    include WithTemplateValidation

    def perform
      return false if @errors.any?

      Magma.instance.db[:models].where(
        project_name: @project_name,
        model_name: @action_params[:model_name],
      ).update(
        template_project_name: clear_template? ? nil : template_project_name,
        template_model_name: clear_template? ? nil : @action_params[:template_model_name],
      )

      true
    end

    private

    def validations
      [
        :validate_model,
        :validate_db_model,
        :validate_template_reference,
        :validate_template_target,
      ]
    end

    def clear_template?
      !!@action_params[:clear_template]
    end

    def validate_template_reference
      return if clear_template?

      return if @action_params[:template_model_name] && @action_params[:template_model_name] != ""

      @errors << Magma::ActionError.new(
        message: "Must include :template_model_name parameter",
        source: @action_params.slice(:action_name, :model_name),
      )
    end

    def validate_template_target
      return if clear_template?
      return unless @errors.empty?

      validate_template_target!(
        self_template_message: "Model cannot point to itself as its template",
        self_template_source_keys: [:action_name, :model_name, :template_project_name, :template_model_name],
        missing_template_message: "Template model does not exist.",
        missing_template_source_keys: [:action_name, :template_project_name, :template_model_name]
      )
    end
  end
end
