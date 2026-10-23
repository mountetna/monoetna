class Magma
  class ProjectToTemplateMapping
    attr_reader :project_name, :template_project

    def initialize(project_name)
      @project_name = project_name
      @template_project = Magma.instance.config(:template_project)
    end

    def models
      @models ||= db[:models].
        where(project_name: @project_name).
        exclude(model_name: 'project').
        select(:project_name, :model_name, :template_project_name, :template_model_name).
        order(:model_name).
        all
    end

    def valid_models
      partition_models[:valid]
    end

    def unmapped_models
      partition_models[:unmapped]
    end

    def invalid_models
      partition_models[:invalid]
    end

    def child_model_name(template_model_name)
      valid_models.find { |model| model[:template_model_name] == template_model_name }&.dig(:model_name)
    end

    def template_model_for(child_model_name)
      valid_models.find { |model| model[:model_name] == child_model_name }&.dig(:template_model_name)
    end

    def project_attributes(model_name)
      @project_attributes ||= db[:attributes].
        where(project_name: @project_name).
        select(:model_name, :attribute_name, :column_name, :type, :link_model_name, :validation).
        all.
        group_by { |attribute| attribute[:model_name] }

      @project_attributes.fetch(model_name, []).
        to_h { |attribute| [attribute[:attribute_name], attribute] }
    end

    def table_name(model_name)
      Sequel[@project_name.to_sym][model_name.pluralize.to_sym]
    end

    def table(model_name)
      db[table_name(model_name)]
    end

    def db
      Magma.instance.db
    end

    private

    def partition_models
      @partition_models ||= begin
        unmapped, mapped = models.partition do |model|
          blank?(model[:template_project_name]) && blank?(model[:template_model_name])
        end
        valid, invalid = mapped.partition do |model|
          model[:template_project_name] == @template_project &&
            template_model_names.include?(model[:template_model_name])
        end

        { valid: valid, unmapped: unmapped, invalid: invalid }
      end
    end

    def template_model_names
      @template_model_names ||= db[:models].
        where(project_name: @template_project).
        select_map(:model_name)
    end

    def blank?(value)
      value.nil? || value.empty?
    end
  end
end
