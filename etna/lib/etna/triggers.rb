module Etna
  module Triggers
    def self.included(host_class)
      host_class.extend ClassMethods
    end

    module ClassMethods
      def before_insert_trigger(table, name, handler_code)
        @triggers ||= []
        @triggers << apply_insert_trigger_sql(:before, table, name, handler_code)
      end

      def before_update_trigger(attrs, table, name, handler_code)
        @triggers ||= []
        @triggers << apply_update_trigger_sql(:before, attrs, table, name, handler_code)
      end

      def before_delete_trigger(table, name, handler_code)
        @triggers ||= []
        @triggers << apply_delete_trigger_sql(:before, table, name, handler_code)
      end

      def triggers_sql
        @triggers.join(';')
      end

      def run_triggers
        self.instance.db.transaction do
          self.instance.db.execute(triggers_sql)
        end
      end
    def apply_insert_trigger_sql(pos, table, name, handler_code)
      apply_trigger_sql(pos, :insert, [], table, name, handler_code)
    end

    def apply_delete_trigger_sql(pos, table, name, handler_code)
      apply_trigger_sql(pos, :delete, [], table, name, handler_code)
    end

    def apply_update_trigger_sql(pos, attrs, table, name, handler_code)
      apply_trigger_sql(pos, :update, attrs, table, name, handler_code)
    end

    def apply_trigger_sql(pos, action, attrs, table, name, handler_code)
      @trigger_sql_cache ||= {}
      attrs_desc = attrs.empty? ? "" : "_when_#{attrs.join('_or_')}_change"
      full_name = "#{pos}_#{table}_#{action}#{attrs_desc}_#{name}"
      if @trigger_sql_cache.include?(full_name)
        unless @trigger_sql_cache[full_name] == handler_code
          raise "Duplication definition for trigger #{full_name}"
        end
      end

      when_sql = attrs.empty? ? "" : "WHEN (#{attrs.map {|a| "(OLD.#{a} IS DISTINCT FROM NEW.#{a})"}.join(' OR ')})"
      @trigger_sql_cache[full_name] ||= <<-SQL
CREATE OR REPLACE FUNCTION handle_#{full_name}() 
RETURNS TRIGGER
  AS $$
BEGIN
   #{handler_code}
END;
$$ LANGUAGE PLPGSQL;

DROP TRIGGER IF EXISTS #{full_name} ON #{table};

CREATE TRIGGER #{full_name}
  #{pos} #{action} ON #{table}
  FOR EACH ROW
#{when_sql}
  EXECUTE PROCEDURE handle_#{full_name}();
SQL
    end
      end
  end
end