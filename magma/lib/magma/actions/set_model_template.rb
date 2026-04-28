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
        :validate_set_template_reference,
      ]
    end

    def clear_template?
      !!@action_params[:clear_template]
    end

    def validate_set_template_reference
      return if clear_template?
      return unless @errors.empty? || !template_model_name_provided?

      validate_template_reference!
    end
  end
end
