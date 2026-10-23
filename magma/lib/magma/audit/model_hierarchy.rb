class Magma
  class ModelHierarchy
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
