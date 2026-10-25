class Magma
  class ModelHierarchy
    # The project model is the root and has no parent. Landing on it is not a subject.
    ROOT_MODEL = 'project'

    attr_reader :project_name

    def initialize(project_name)
      @project_name = project_name
    end

    # { model_name => parent_model_name }
    def parents
      @parents ||= attributes_of_type('parent').to_h do |attribute|
        [attribute[:model_name], attribute[:link_model_name]]
      end
    end

    def parent_attribute(model_name)
      attributes_of_type('parent').find { |attribute| attribute[:model_name] == model_name }
    end

    # Parent steps from start_model up to, but not including, target_model.
    # Each step is { model:, column: } for that model's parent attribute.
    # Nil when the walk cycles, loses its parent, hits the project root, or never reaches target.
    def parent_path(start_model, target_model)
      return if start_model.nil? || target_model.nil? || target_model == ROOT_MODEL

      steps = []
      current = start_model
      seen = {}

      until current == target_model
        return if current.nil? || current == ROOT_MODEL || seen[current]

        seen[current] = true
        attribute = parent_attribute(current)
        parent_name = attribute && attribute[:link_model_name]
        return if parent_name.nil?

        steps << { model: current, column: attribute[:column_name] }
        current = parent_name
      end

      steps
    end

    def link_between(model_name, target_model)
      attributes_of_type('link').find do |attribute|
        attribute[:model_name] == model_name && attribute[:link_model_name] == target_model
      end
    end

    private

    def attributes_of_type(type)
      attributes.select { |attribute| attribute[:type] == type }
    end

    def attributes
      @attributes ||= Magma.instance.db[:attributes].
        where(project_name: @project_name).
        select(:model_name, :attribute_name, :column_name, :type, :link_model_name).
        all
    end
  end
end
