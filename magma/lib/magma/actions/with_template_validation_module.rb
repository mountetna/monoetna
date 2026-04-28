class Magma
  module WithTemplateValidation
    DEFAULT_TEMPLATE_PROJECT = "coprojects_template"

    def validate_template_reference!
      unless template_model_name_provided?
        add_missing_template_model_name_error
        return
      end

      unless template_not_self_reference?
        add_template_self_reference_error
        return
      end

      return if referenced_template_model_exists?

      add_missing_template_target_error
    end

    private

    def template_project_name
      @action_params[:template_project_name] || DEFAULT_TEMPLATE_PROJECT
    end

    def template_model_name_provided?
      @action_params[:template_model_name] && @action_params[:template_model_name] != ""
    end

    def template_not_self_reference?
      @project_name.to_s != template_project_name.to_s ||
        @action_params[:model_name].to_s != @action_params[:template_model_name].to_s
    end

    def referenced_template_model_exists?
      Magma.instance.db[:models].where(
        project_name: template_project_name,
        model_name: @action_params[:template_model_name]
      ).first
    end

    def add_missing_template_model_name_error
      @errors << Magma::ActionError.new(
        message: "Must include :template_model_name parameter",
        source: @action_params.slice(:action_name, :model_name),
      )
    end

    def add_template_self_reference_error
      @errors << Magma::ActionError.new(
        message: "Model cannot point to itself as its template",
        source: @action_params.slice(:action_name, :model_name, :template_project_name, :template_model_name),
      )
    end

    def add_missing_template_target_error
      @errors << Magma::ActionError.new(
        message: "Template model does not exist.",
        source: @action_params.slice(:action_name, :template_project_name, :template_model_name),
      )
    end
  end
end
