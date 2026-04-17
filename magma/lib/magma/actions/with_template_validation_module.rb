class Magma
  module WithTemplateValidation
    DEFAULT_TEMPLATE_PROJECT = "coprojects_template"

    private

    def template_project_name
      @action_params[:template_project_name] || DEFAULT_TEMPLATE_PROJECT
    end

    def validate_template_target!(
      self_template_message:,
      self_template_source_keys:,
      missing_template_message:,
      missing_template_source_keys:
    )
      if template_self_reference?
        @errors << Magma::ActionError.new(
          message: self_template_message,
          source: @action_params.slice(*self_template_source_keys),
        )
        return
      end

      return if template_model_exists?

      @errors << Magma::ActionError.new(
        message: missing_template_message,
        source: @action_params.slice(*missing_template_source_keys),
      )
    end

    def template_self_reference?
      @project_name.to_s == template_project_name.to_s &&
        @action_params[:model_name].to_s == @action_params[:template_model_name].to_s
    end

    def template_model_exists?
      Magma.instance.db[:models].where(
        project_name: template_project_name,
        model_name: @action_params[:template_model_name]
      ).first
    end
  end
end
