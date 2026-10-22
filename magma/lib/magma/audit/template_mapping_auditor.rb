class Magma
  class TemplateMappingAuditor
    def initialize(mapping)
      @mapping = mapping
      @template_project = mapping.template_project
    end

    def report
      project_name = @mapping.project_name
      unmapped, valid, invalid = @mapping.unmapped_models, @mapping.valid_models, @mapping.invalid_models

      issues = {
        unmapped_models: unmapped.map { |model| model[:model_name] },
        invalid_mappings: invalid.map { |model| model_mapping(model) },
        missing_template_columns: valid.filter_map { |model| check_required_on_child_models(project_name, model) },
        ontology_mismatches: valid.flat_map { |model| check_ontology_mismatches(project_name, model) },
        invalid_ontology_values: valid.flat_map { |model| check_ontology_values(project_name, model) }
      }

      {
        conforming: issues.values.all?(&:empty?),
        summary: audit_summary(project_name, issues),
        mapped_models: valid.map { |model| model_mapping(model) },
        **issues,
        required_column_results: valid.flat_map { |model| required_column_results(project_name, model) },
        ontology_value_results: valid.flat_map { |model| ontology_value_results(project_name, model) },
        missing_ontology_columns: valid.flat_map { |model| missing_ontology_columns(project_name, model) }
      }
    end

    private

    def audit_summary(project_name, issues)
      valid = @mapping.valid_models
      required_total = valid.sum { |model| required_attributes(model).length }
      required_missing = issues[:missing_template_columns].sum { |entry| entry[:columns].length }
      ontology_total = valid.sum { |model| required_ontology_attributes(model).length }
      ontology_missing = valid.sum do |model|
        columns = child_columns(project_name, model)
        required_ontology_attributes(model).count { |attribute| !columns.key?(attribute[:attribute_name]) }
      end

      {
        models: {
          total: @mapping.models.length,
          mapped: valid.length,
          unmapped: @mapping.unmapped_models.length,
          invalid: @mapping.invalid_models.length
        },
        required_columns: {
          total: required_total,
          confirmed: required_total - required_missing,
          missing: required_missing
        },
        ontologies: {
          total: ontology_total,
          confirmed: ontology_total - ontology_missing - issues[:ontology_mismatches].length,
          missing_columns: ontology_missing,
          mismatched: issues[:ontology_mismatches].length,
          invalid_values: issues[:invalid_ontology_values].sum { |entry| entry[:values].length }
        }
      }
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

    def required_column_results(project_name, model)
      columns = child_columns(project_name, model)

      required_attributes(model).map do |attribute|
        model_mapping(model).merge(
          column: attribute[:attribute_name],
          present: columns.key?(attribute[:attribute_name])
        )
      end
    end

    def ontology_value_results(project_name, model)
      required_ontology_columns(project_name, model).flat_map do |attribute, ontology, column|
        next [] unless ontology_validation(column) == ontology

        validator = ontology_validator(ontology)
        child_column_values(project_name, model, column[:column_name]).map do |value|
          model_mapping(model).merge(
            column: attribute[:attribute_name],
            value: value,
            valid: validator.validate(value)
          )
        end
      end
    end

    def missing_ontology_columns(project_name, model)
      columns = child_columns(project_name, model)

      required_ontology_attributes(model).filter_map do |attribute|
        next if columns.key?(attribute[:attribute_name])

        model_mapping(model).merge(column: attribute[:attribute_name])
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

    def required_ontology_attributes(model)
      required_attributes(model).select { |attribute| ontology_validation(attribute) }
    end

    def ontology_validation(attribute)
      validation = attribute[:validation]&.to_h&.symbolize_keys
      validation if validation && validation[:type] == 'Ontology'
    end

    def ontology_validator(ontology)
      @ontology_validators ||= {}
      @ontology_validators[ontology] ||= Magma::ValidationObject.build(ontology)
    end

    def model_mapping(model)
      { model: model[:model_name], template_model: model[:template_model_name] }
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
      @mapping.project_attributes(model[:model_name])
    end

    def child_column_values(project_name, model, column)
      Magma.instance.db[Sequel[project_name.to_sym][model[:model_name].pluralize.to_sym]].
        exclude(column.to_sym => nil).
        exclude(column.to_sym => '').
        distinct.
        select_order_map(column.to_sym)
    end
  end
end
