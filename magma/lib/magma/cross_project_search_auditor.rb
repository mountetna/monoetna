class Magma
  class CrossProjectSearchAuditor
    def initialize(project_name: nil)
      @project_name = project_name
      @template_project = Magma.instance.config(:template_project)
    end

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
        model[:template_project_name] == @template_project &&
          template_model_names.include?(model[:template_model_name])
      end

      issues = {
        unmapped_models: unmapped.map { |model| model[:model_name] },
        invalid_mappings: invalid.map { |model| model_mapping(model) },
        missing_template_columns: valid.filter_map { |model| check_required_on_child_models(project_name, model) },
        ontology_mismatches: valid.flat_map { |model| check_ontology_mismatches(project_name, model) },
        invalid_ontology_values: valid.flat_map { |model| check_ontology_values(project_name, model) }
      }

      { project: project_name, conforming: issues.values.all?(&:empty?), **issues }
    end

    def check_required_on_child_models(project_name, model)
      missing = required_attributes(model).map { |attribute| attribute[:attribute_name] } - child_columns(project_name, model).keys
      return if missing.empty?

      model_mapping(model).merge(columns: missing)
    end

    def check_ontology_mismatches(project_name, model)
      required_ontology_columns(project_name, model).filter_map do |attribute, ontology, column|
        actual = ontology_validation(column)
        next if actual == ontology

        model_mapping(model).merge(
          column: attribute[:attribute_name],
          expected: ontology.except(:type),
          actual: actual&.except(:type)
        )
      end
    end

    def check_ontology_values(project_name, model)
      required_ontology_columns(project_name, model).filter_map do |attribute, ontology, column|
        validator = ontology_validator(ontology)
        values = child_column_values(project_name, model, column[:column_name]).reject { |value| validator.validate(value) }
        next if values.empty?

        model_mapping(model).merge(column: attribute[:attribute_name], ontology: ontology.except(:type), values: values)
      end
    end

    def required_attributes(model)
      template_attributes(model).select { |attribute| attribute[:required_on_child_models] }
    end

    def required_ontology_columns(project_name, model)
      columns = child_columns(project_name, model)

      required_attributes(model).filter_map do |attribute|
        ontology = ontology_validation(attribute)
        column = columns[attribute[:attribute_name]]
        [attribute, ontology, column] if ontology && column
      end
    end

    def ontology_validation(attribute)
      validation = attribute[:validation]&.to_h&.symbolize_keys
      validation if validation && validation[:type] == 'Ontology'
    end

    def ontology_validator(ontology)
      @ontology_validators ||= {}
      @ontology_validators[ontology] ||= Magma::ValidationObject.build(ontology)
    end

    def excluded_projects
      @excluded_projects ||= [@template_project] + all_template_attributes.filter_map do |attribute|
        ontology_validation(attribute)&.dig(:project_name)
      end.uniq
    end

    def model_mapping(model)
      { model: model[:model_name], template_model: model[:template_model_name] }
    end

    def project_models
      models = Magma.instance.db[:models].
        exclude(project_name: excluded_projects).
        exclude(model_name: 'project')
      models = models.where(project_name: @project_name) if @project_name

      models.
        select(:project_name, :model_name, :template_project_name, :template_model_name).
        order(:project_name, :model_name).
        all
    end

    def template_model_names
      @template_model_names ||= Magma.instance.db[:models].
        where(project_name: @template_project).
        select_map(:model_name)
    end

    def template_attributes(model)
      @template_attributes ||= all_template_attributes.group_by { |attribute| attribute[:model_name] }
      @template_attributes.fetch(model[:template_model_name], [])
    end

    def all_template_attributes
      @all_template_attributes ||= Magma.instance.db[:attributes].
        where(project_name: @template_project).
        select(:model_name, :attribute_name, :required_on_child_models, :validation).
        order(:attribute_name).
        all
    end

    def child_columns(project_name, model)
      @child_columns ||= Magma.instance.db[:attributes].
        exclude(project_name: excluded_projects).
        select(:project_name, :model_name, :attribute_name, :column_name, :validation).
        all.
        group_by { |attribute| [attribute[:project_name], attribute[:model_name]] }

      @child_columns.fetch([project_name, model[:model_name]], []).
        to_h { |attribute| [attribute[:attribute_name], attribute] }
    end

    def child_column_values(project_name, model, column)
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
