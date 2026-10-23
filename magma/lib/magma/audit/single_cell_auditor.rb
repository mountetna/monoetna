require_relative 'model_hierarchy'

class Magma
  class SingleCellAuditor
    # Root to leaf. sc_seq also links to sc_seq_pool, which is not a parent.
    SINGLE_CELL_HIERARCHY = %w[subject biospecimen sc_seq].freeze
    SINGLE_CELL_POOL = 'sc_seq_pool'

    def initialize(mapping)
      @mapping = mapping
      @child_hierarchy = ModelHierarchy.new(mapping.project_name)
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
      parent_column = subject_path&.first&.dig(:column)&.to_sym
      depooled = linked && parent_column && linked.exclude(parent_column => nil)
      # Parent column null, with or without a pool. pool_linked - depooled misses rows that have neither.
      missing_subject = parent_column && sc_seq_table.where(parent_column => nil).count

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
            depooled: depooled&.count,
            missing_subject: missing_subject
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
