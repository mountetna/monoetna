Sequel::Schema::AlterTableGenerator.class_eval do
  def drop_foreign_constraint_if_exists(column)
    @operations << {:op => :drop_if_exists_constraint, :type => :foreign_key, :columns => [column]}
    nil
  end
end

Sequel::Database.class_eval do
  def alter_table_drop_if_exists_constraint_sql(table, op)
    orig = alter_table_drop_constraint_sql(table, op)
    l, r = orig.split("CONSTRAINT")
    "#{l}CONSTRAINT IF EXISTS#{r}"
  rescue Sequel::Error => e
    # If e is a "missing foreign key" then there is no way to construct
    #   the constraint name, because a foreign key does not exist.
    # So skip the drop_constraint in this case.
    raise e unless e.message =~ /^Missing foreign key for/
    nil
  end
end

class Magma
  class Migration
    class << self
      def table_name(model)
        "Sequel[:#{model.project_name}][:#{model.implicit_table_name}]"
      end
      def create(model)
        if Magma.instance.db.table_exists?(model.table_name)
          return Magma::UpdateMigration.new(model)
        else
          return Magma::CreateMigration.new(model)
        end
      end
    end

    def initialize(model)
      @model = model
      @changes = {}
    end

    def change(key, lines)
      @changes[key] ||= []
      @changes[key].concat lines
    end

    def empty?
      @changes.empty?
    end

    def to_s
      @changes.map do |key,lines|
        lines.empty? ?
          space(key, 2) :
          [
            space("#{key} do", 2),
            lines.map do |line|
              space(line,3)
            end,
            space('end',2)
          ].flatten.join("\n")
      end.join("\n").chomp
    end

    private

    def attribute_migration(att)
      case att
      when Magma::ForeignKeyAttribute
        [
          foreign_key_entry(att.column_name, att.link_model),
          index_entry(att.column_name)
        ]
      when Magma::Attribute
        [
          column_entry(att.column_name, att.database_type),
          att.unique && unique_entry(att.column_name),
          att.index && index_entry(att.column_name)
        ].compact
      else
        nil
      end
    end

    # this denotes a link attribute that points here (i.e.,
    # the foreign key is in the link_model) and therefore has
    # no column in this model
    def foreign_attribute?(att)
      [ Magma::ChildAttribute, Magma::CollectionAttribute, Magma::TableAttribute ].any? do |att_class|
        att.instance_of?(att_class)
      end
    end

    def schema_supports_attribute?(model, att)
      if foreign_attribute?(att)
        return true
      else
        return model.schema.has_key?(att.column_name.to_sym)
      end
    end

    def schema_unchanged?(model, att)
      # we don't need to worry about models that link to us
      return true if foreign_attribute?(att)
      # neither can foreign keys change their type
      return true if att.is_a?(Magma::ForeignKeyAttribute)

      literal_type = att.is_a?(DateTimeAttribute)?  :"timestamp without time zone" :
        Magma.instance.db.cast_type_literal(att.database_type)

      return model.schema[att.column_name.to_sym][:db_type].to_sym == literal_type.to_sym
    end

    SPC='  '
    def space(txt, pad)
      "#{SPC*pad}#{txt}"
    end

    def now
      DateTime.now.to_time.to_i
    end
  end
  class CreateMigration < Migration
    def initialize(model)
      super
      tlb_nm = "create_table(#{Magma::Migration.table_name(model)})"
      change(tlb_nm, ['primary_key :id']+new_attributes)
    end

    def new_attributes
      @model.attributes.reject { |name, attr| attr.primary_key? }.map do |name,att|
        next if foreign_attribute?(att)
        attribute_migration(att)
      end.compact.flatten
    end

    def foreign_key_entry column_name, foreign_model
      "foreign_key :#{column_name}, #{Magma::Migration.table_name(foreign_model)}"
    end

    def column_entry name, type
      "#{type.respond_to?(:name) ? type.name : type} :#{name}"
    end

    def unique_entry name
      "unique :#{name}"
    end

    def index_entry column
      if column.is_a? Array
        "index [#{column.map{|c| ":#{c}"}.join(", ")}]"
      else
        "index :#{column}"
      end
    end
  end

  class UpdateMigration < Migration
    def initialize model
      super
      change("alter_table(#{Magma::Migration.table_name(model)})", missing_attributes) unless missing_attributes.empty?

      change("alter_table(#{Magma::Migration.table_name(model)})", removed_attributes) unless removed_attributes.empty?

      change("alter_table(#{Magma::Migration.table_name(model)})", changed_attributes) unless changed_attributes.empty?
    end

    def foreign_key_entry column_name, foreign_model
      "add_foreign_key :#{column_name}, #{Magma::Migration.table_name(foreign_model)}"
    end

    def column_type_entry name, type
      "set_column_type :#{name}, #{type.is_a?(Symbol) ? ":#{type}, using: '#{name}::#{type}'" : type}"
    end

    def column_entry name, type
      "add_column :#{name}, #{type.is_a?(Symbol) ? ":#{type}" : type}"
    end

    def remove_column_entry name
      new_name = "#{name}_#{now}_backup"

      # Drop constraints so can add new records without specifying
      #   a FK.
      [
        "drop_index :#{name}, if_exists: true",
        "drop_foreign_constraint_if_exists :#{name}",
        "rename_column :#{name}, :#{new_name}"
      ]
    end

    def unique_entry name
      "add_unique_constraint :#{name}"
    end

    def index_entry column
      if column.is_a? Array
        "add_index [#{column.map{|c| ":#{c}"}.join(", ")}]"
      else
        "add_index :#{column}"
      end
    end

    private

    def missing_attributes
      @missing_attributes ||= @model.attributes.reject { |name, attr| attr.primary_key? }.map do |name,att|
        next if schema_supports_attribute?(@model, att)
        attribute_migration(att)
      end.compact.flatten
    end


    def changed_attributes
      @changed_attributes ||= @model.attributes.reject { |name, attr| attr.primary_key? }.map do |name,att|
        next unless schema_supports_attribute?(@model,att)
        next if schema_unchanged?(@model,att)
        column_type_entry(att.column_name, att.database_type)
      end.compact.flatten
    end

    def removed_attributes
      @removed_attributes ||= @model.schema.map do |name, db_opts|
        next if @model.attributes[name]
        next if @model.attributes[ name.to_s.sub(/_id$/,'').to_sym ]
        next if @model.attributes[ name.to_s.sub(/_type$/,'').to_sym ]
        next if name =~ /\d+_backup$/
        next if db_opts[:primary_key]
        next if attribute_with_different_column_name?(name)
        remove_column_entry(name)
      end.compact.flatten
    end

    def attribute_with_different_column_name?(name)
      @model.attributes.any? do |attribute_name, attribute|
        attribute.column_name.to_sym == name
      end
    end
  end

  class RemoveModelMigration < Migration
    def initialize(model)
      super

      drop_index = "alter_table(#{Magma::Migration.table_name(model)})"
      change(drop_index, remove_indices)

      new_name = "Sequel[:#{model.project_name}][:#{model.implicit_table_name}_#{now}_backup]"
      table_rename = "rename_table(#{Magma::Migration.table_name(model)}, #{new_name})"
      change(table_rename, [])
    end

    private

    def remove_indices
      # Need to remove the FK index so can
      #   re-add model with the same name, if needed
      [
        "drop_index [:#{parent_attribute.column_name}]"
      ]
    end

    def parent_attribute
      @model.attributes.values.select do |attribute|
        attribute.is_a?(Magma::ParentAttribute)
      end.first
    end
  end
end
