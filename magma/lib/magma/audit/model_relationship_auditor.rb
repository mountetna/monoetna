require_relative 'model_hierarchy'
require_relative 'single_cell_auditor'

class Magma
  class ModelRelationshipAuditor
    def initialize(mapping)
      @mapping = mapping
      @child_hierarchy = ModelHierarchy.new(mapping.project_name)
      @template_hierarchy = ModelHierarchy.new(mapping.template_project)
    end

    def hierarchy
      expected = expected_parent_edges
      actual = actual_parent_edges

      issues = (expected.keys | actual.keys).sort.filter_map do |model|
        next if expected[model] == actual[model]

        {
          model: model,
          template_model: @mapping.template_model_for(model),
          parent: actual[model],
          expected_parent: expected[model]
        }
      end

      { conforming: issues.empty?, issues: issues }
    end

    def single_cell
      SingleCellAuditor.new(@mapping).single_cell
    end

    private

    def expected_parent_edges
      @template_hierarchy.parents.each_with_object({}) do |(template_model, template_parent), edges|
        model = @mapping.child_model_name(template_model)
        parent = @mapping.child_model_name(template_parent)
        edges[model] = parent if model && parent
      end
    end

    def actual_parent_edges
      @child_hierarchy.parents.each_with_object({}) do |(model, parent), edges|
        next unless @mapping.template_model_for(model) && @mapping.template_model_for(parent)

        edges[model] = parent
      end
    end
  end
end
