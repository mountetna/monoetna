require_relative 'model_hierarchy'

class Magma
  class ModelRelationshipAuditor
    # Root to leaf. sc_seq also links to sc_seq_pool, which is not a parent.
    SINGLE_CELL_HIERARCHY = %w[subject biospecimen sc_seq].freeze
    SINGLE_CELL_POOL = 'sc_seq_pool'

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
      child_models = single_cell_child_models
      missing = child_models.select { |_, name| name.nil? }.keys

      unless missing.empty?
        return {
          applicable: false,
          conforming: true,
          reason: "Child project has no model for #{missing.join(', ')}",
          child_models: child_models
        }
      end

      sc_seq = child_models[:sc_seq]
      pool = child_models[:sc_seq_pool]
      pool_link = @child_hierarchy.link_between(sc_seq, pool)
      subject_path = single_cell_subject_path(child_models)

      issues = []
      issues << "#{sc_seq} has no link attribute to #{pool}" unless pool_link
      issues << "#{sc_seq} parent chain does not follow #{SINGLE_CELL_HIERARCHY.join(' -> ')}" unless subject_path

      sc_seq_table = @mapping.table(sc_seq)
      linked = pool_link && sc_seq_table.exclude(pool_link[:column_name].to_sym => nil)
      depooled = linked && subject_path && linked.exclude(subject_path.first[:column].to_sym => nil)

      {
        applicable: true,
        conforming: issues.empty?,
        child_models: child_models,
        pool_link: pool_link&.dig(:attribute_name),
        subject_path: subject_path&.map { |step| step[:model] }&.push(child_models[:subject]),
        issues: issues,
        counts: {
          pools: {
            total: @mapping.table(pool).count,
            depooled: linked&.select(pool_link[:column_name].to_sym)&.distinct&.count
          },
          sc_seq: {
            total: sc_seq_table.count,
            pool_linked: linked&.count,
            depooled: depooled&.count
          },
          subjects: depooled && depooled_subject_count(pool_link, subject_path)
        }
      }
    end

    private

    def single_cell_child_models
      (SINGLE_CELL_HIERARCHY + [SINGLE_CELL_POOL]).to_h do |template_model|
        [template_model.to_sym, @mapping.child_model_name(template_model)]
      end
    end

    # Leaf to root, so the first step is sc_seq. Nil unless every parent matches.
    def single_cell_subject_path(child_models)
      path = []
      SINGLE_CELL_HIERARCHY.reverse.each_cons(2) do |template_model, template_parent|
        model = child_models[template_model.to_sym]
        parent = child_models[template_parent.to_sym]
        attribute = @child_hierarchy.parent_attribute(model)
        return if attribute.nil? || attribute[:link_model_name] != parent

        path << { model: model, column: attribute[:column_name] }
      end
      path
    end

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

    def depooled_subject_count(pool_link, subject_path)
      aliases = subject_path.each_index.map { |i| :"t#{i}" }
      dataset = @mapping.db.from(Sequel.as(@mapping.table_name(subject_path.first[:model]), aliases.first)).
        exclude(Sequel[aliases.first][pool_link[:column_name].to_sym] => nil)

      subject_path.each_cons(2).with_index do |(child, parent), i|
        dataset = dataset.join(
          Sequel.as(@mapping.table_name(parent[:model]), aliases[i + 1]),
          id: Sequel[aliases[i]][child[:column].to_sym]
        )
      end

      subject_column = Sequel[aliases.last][subject_path.last[:column].to_sym]
      dataset.exclude(subject_column => nil).select(subject_column).distinct.count
    end
  end
end
