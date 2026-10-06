class Magma
  class TemplateAudit
    TEMPLATE_PROJECT = 'coprojects_template'.freeze
    EXCLUDED_PROJECTS = [TEMPLATE_PROJECT, Magma::OntologyValidationObject::PROJECT].freeze

    def report
      projects = project_models.group_by { |model| model[:project_name] }

      { projects: projects.map { |project_name, models| audit_project(project_name, models) } }
    end

    private

    def audit_project(project_name, models)
      unmapped, mapped = models.partition do |model|
        blank?(model[:template_project_name]) && blank?(model[:template_model_name])
      end
      valid, invalid = mapped.partition do |model|
        model[:template_project_name] == TEMPLATE_PROJECT &&
          template_model_names.include?(model[:template_model_name])
      end

      issues = {
        unmapped_models: unmapped.map { |model| model[:model_name] },
        invalid_mappings: invalid.map { |model| model_mapping(model) },
        missing_template_columns: valid.filter_map { |model| missing_columns(project_name, model) },
        invalid_ontology_values: valid.flat_map { |model| invalid_values(project_name, model) }
      }

      { project: project_name, conforming: issues.values.all?(&:empty?), **issues }
    end

    def missing_columns(project_name, model)
      required = template_attributes(model).select { |attribute| attribute[:template_required] }
      missing = required.map { |attribute| attribute[:attribute_name] } - local_columns(project_name, model).keys
      return if missing.empty?

      model_mapping(model).merge(columns: missing)
    end

    def invalid_values(project_name, model)
      template_attributes(model).filter_map do |attribute|
        next unless attribute[:validation_type] == 'Ontology'

        column = local_columns(project_name, model)[attribute[:attribute_name]]
        next unless column

        table = attribute[:validation_value]
        validator = Magma::OntologyValidationObject.new(value: table)
        values = stored_values(project_name, model, column).reject { |value| validator.validate(value) }
        next if values.empty?

        model_mapping(model).merge(column: attribute[:attribute_name], table: table, values: values)
      end
    end

    def model_mapping(model)
      { model: model[:model_name], template_model: model[:template_model_name] }
    end

    def project_models
      Magma.instance.db[:models].
        exclude(project_name: EXCLUDED_PROJECTS).
        select(:project_name, :model_name, :template_project_name, :template_model_name).
        order(:project_name, :model_name).
        all
    end

    def template_model_names
      @template_model_names ||= Magma.instance.db[:models].
        where(project_name: TEMPLATE_PROJECT).
        select_map(:model_name)
    end

    def template_attributes(model)
      validation = Sequel.pg_json_op(:validation)

      @template_attributes ||= Magma.instance.db[:attributes].
        where(project_name: TEMPLATE_PROJECT).
        select(
          :model_name,
          :attribute_name,
          :template_required,
          validation.get_text('type').as(:validation_type),
          validation.get_text('value').as(:validation_value)
        ).
        order(:attribute_name).
        all.
        group_by { |attribute| attribute[:model_name] }

      @template_attributes.fetch(model[:template_model_name], [])
    end

    def local_columns(project_name, model)
      @local_columns ||= Magma.instance.db[:attributes].
        exclude(project_name: EXCLUDED_PROJECTS).
        select(:project_name, :model_name, :attribute_name, :column_name).
        all.
        group_by { |attribute| [attribute[:project_name], attribute[:model_name]] }

      @local_columns.fetch([project_name, model[:model_name]], []).
        to_h { |attribute| [attribute[:attribute_name], attribute[:column_name]] }
    end

    def stored_values(project_name, model, column)
      Magma.instance.db[Sequel[project_name.to_sym][model[:model_name].pluralize.to_sym]].
        exclude(column.to_sym => nil).
        exclude(column.to_sym => '').
        distinct.
        select_order_map(column.to_sym)
    end

    def blank?(value)
      value.nil? || value.empty?
    end
  end
end
