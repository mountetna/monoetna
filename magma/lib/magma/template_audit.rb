class Magma
  class TemplateAudit
    TEMPLATE_PROJECT = 'coprojects_template'.freeze

    def call
      { projects: audited_projects }
    end

    private

    def audited_projects
      models_by_project.keys.sort.map do |project_name|
        audit_project(project_name, models_by_project[project_name])
      end
    end

    def models_by_project
      @models_by_project ||= project_models.group_by { |model| model[:project_name] }
    end

    def audit_project(project_name, models)
      model_groups = classify_models(models)

      issues = {
        unmapped_models: unmapped_model_names(model_groups[:unmapped]),
        invalid_mappings: invalid_mapping_reports(model_groups[:invalid]),
        missing_template_columns: missing_template_columns(project_name, model_groups[:valid])
      }

      {
        project: project_name,
        conforming: issues.values.all?(&:empty?),
        **issues
      }
    end

    def classify_models(models)
      unmapped, mapped = models.partition { |model| mapping_blank?(model) }
      valid, invalid = mapped.partition { |model| valid_mapping?(model) }

      {
        unmapped: unmapped,
        valid: valid,
        invalid: invalid
      }
    end

    def unmapped_model_names(models)
      models.map { |model| model[:model_name] }.sort
    end

    def invalid_mapping_reports(models)
      models.map do |model|
        {
          model: model[:model_name],
          template_model: model[:template_model_name]
        }
      end.sort_by { |mapping| mapping[:model] }
    end

    def missing_template_columns(project_name, models)
      models.filter_map do |model|
        expected = enforced_attributes_for(model[:template_model_name])
        actual = attributes_for(project_name, model[:model_name])
        missing = expected.reject { |attribute_name| actual.include?(attribute_name) }
        next if missing.empty?

        {
          model: model[:model_name],
          template_model: model[:template_model_name],
          columns: missing
        }
      end
    end

    def mapping_blank?(model)
      blank?(model[:template_project_name]) && blank?(model[:template_model_name])
    end

    def valid_mapping?(model)
      model[:template_project_name] == TEMPLATE_PROJECT &&
        template_model_names.include?(model[:template_model_name])
    end

    def enforced_attributes_for(model_name)
      template_enforced_attributes.fetch(model_name, [])
    end

    def attributes_for(project_name, model_name)
      local_attributes.fetch([project_name, model_name], [])
    end

    def project_models
      Magma.instance.db[:models].
        exclude(project_name: TEMPLATE_PROJECT).
        select(:project_name, :model_name, :template_project_name, :template_model_name).
        order(:project_name, :model_name).
        all
    end

    def template_model_names
      @template_model_names ||= Magma.instance.db[:models].
        where(project_name: TEMPLATE_PROJECT).
        select_map(:model_name)
    end

    def template_enforced_attributes
      @template_enforced_attributes ||= Magma.instance.db[:attributes].
        where(project_name: TEMPLATE_PROJECT, template_enforced: true).
        select(:model_name, :attribute_name).
        order(:model_name, :attribute_name).
        all.
        group_by { |attribute| attribute[:model_name] }.
        transform_values do |attributes|
          attributes.map { |attribute| attribute[:attribute_name] }
        end
    end

    def local_attributes
      @local_attributes ||= Magma.instance.db[:attributes].
        exclude(project_name: TEMPLATE_PROJECT).
        select(:project_name, :model_name, :attribute_name).
        all.
        group_by { |attribute| [attribute[:project_name], attribute[:model_name]] }.
        transform_values do |attributes|
          attributes.map { |attribute| attribute[:attribute_name] }
        end
    end

    def blank?(value)
      value.nil? || value.empty?
    end
  end
end
