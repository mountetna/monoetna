require_relative 'model_hierarchy'

class Magma
  class SingleCellAuditor
    # Root to leaf. sc_seq also links to sc_seq_pool, which is not a parent.
    # sc_seq must parent to biospecimen. biospecimen may reach subject through other models.
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
      missing_subject = subject_path && missing_subject_count(subject_path)

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

    # Leaf to root, so the first step is sc_seq. Nil unless sc_seq parents to biospecimen
    # and walking parent attributes from biospecimen reaches subject.
    def single_cell_subject_path(child_models)
      sc_seq = child_models[:sc_seq]
      biospecimen = child_models[:biospecimen]
      subject = child_models[:subject]
      sc_seq_parent = @child_hierarchy.parent_attribute(sc_seq)
      return if sc_seq_parent.nil? || sc_seq_parent[:link_model_name] != biospecimen

      upward = @child_hierarchy.parent_path(biospecimen, subject)
      return if upward.nil? || upward.empty?

      [{ model: sc_seq, column: sc_seq_parent[:column_name] }] + upward
    end

    def depooled_subject_count(pool_link, subject_path)
      dataset, aliases = lineage_dataset(subject_path)
      dataset = dataset.exclude(Sequel[aliases.first][pool_link[:column_name].to_sym] => nil)
      subject_column = Sequel[aliases.last][subject_path.last[:column].to_sym]
      dataset.exclude(subject_column => nil).select(subject_column).distinct.count
    end

    # Rows whose parent id is null somewhere between biospecimen and subject.
    # A null biospecimen on sc_seq counts. A dangling id that points at no row does not.
    def missing_subject_count(subject_path)
      dataset, aliases = lineage_dataset(subject_path, left: true)
      conditions = subject_path.each_with_index.map do |step, index|
        column_null = { Sequel[aliases[index]][step[:column].to_sym] => nil }
        next column_null if index.zero?

        Sequel.&(Sequel.negate(Sequel[aliases[index]][:id] => nil), column_null)
      end

      dataset.where(Sequel.|(*conditions)).count
    end

    def lineage_dataset(subject_path, left: false)
      aliases = subject_path.each_index.map { |i| :"t#{i}" }
      dataset = @mapping.db.from(Sequel.as(@mapping.table_name(subject_path.first[:model]), aliases.first))

      subject_path.each_cons(2).with_index do |(child, parent), index|
        join_args = [
          Sequel.as(@mapping.table_name(parent[:model]), aliases[index + 1]),
          { id: Sequel[aliases[index]][child[:column].to_sym] }
        ]
        dataset = left ? dataset.left_join(*join_args) : dataset.join(*join_args)
      end

      [dataset, aliases]
    end
  end
end
