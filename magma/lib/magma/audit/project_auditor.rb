require_relative 'project_to_template_mapping'
require_relative 'template_mapping_auditor'
require_relative 'model_relationship_auditor'

class Magma
  class ProjectAuditor
    def initialize(project_name)
      @mapping = ProjectToTemplateMapping.new(project_name)
    end

    def report
      mapping_auditor = TemplateMappingAuditor.new(@mapping)
      relationships = ModelRelationshipAuditor.new(@mapping)
      model_relationships = {
        hierarchy: relationships.hierarchy,
        single_cell: relationships.single_cell
      }

      {
        project: @mapping.project_name,
        conforming: mapping_auditor.conforming? && model_relationships.values.all? { |check| check[:conforming] },
        template_mapping: mapping_auditor.report,
        model_relationships: model_relationships
      }
    end
  end
end
