require_relative 'project_to_template_mapping'
require_relative 'template_mapping_auditor'
require_relative 'model_relationship_auditor'

class Magma
  class ProjectAuditor
    def initialize(project_name)
      @mapping = ProjectToTemplateMapping.new(project_name)
    end

    def report
      template_mapping = TemplateMappingAuditor.new(@mapping).report
      relationships = ModelRelationshipAuditor.new(@mapping)
      model_relationships = {
        hierarchy: relationships.hierarchy,
        single_cell: relationships.single_cell
      }

      {
        project: @mapping.project_name,
        conforming: template_mapping[:conforming] && model_relationships.values.all? { |check| check[:conforming] },
        template_mapping: template_mapping,
        model_relationships: model_relationships
      }
    end
  end
end
