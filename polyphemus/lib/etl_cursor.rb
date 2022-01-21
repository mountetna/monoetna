require 'date'
require 'yaml'

class Polyphemus
  # An EtlCursor represents a named "Cursor", generally encoding things such as an updated_at, limit, and offset,
  # which can correctly capture the entire state of query for batched etl jobs.
  class EtlCursor
    attr_reader :value, :name
    attr_accessor :updated_at

    def initialize(name, updated_at = Time.at(0), value = {}, version = 0)
      @name = name
      @value = value
      @updated_at = updated_at
      @version = version
    end

    def to_s
      value.inspect
    end

    def [](k)
      value[k.to_s]
    end

    def []=(k, v)
      value[k.to_s] = v
    end

    def self.from_env(env = ENV, prefix='ETL_CURSOR_')
      value = self.class.hash_from_env(env, prefix)
      new_k_params = {}
      new_p_params = []

      self.class.method(:new).parameters.each do |type, p_key|
        if type == :req || type == :opt
          if value.include?(p_key.to_s)
            new_p_params << value[p_key.to_s]
          elsif type == :req
            raise "Value #{p_key} must be set in the environment as #{prefix}#{p_key.to_s.upcase}"
          end
        elsif type == :keyreq || type == :key
          if value.include?(p_key.to_s)
            new_k_params[p_key] = value[p_key.to_s]
          elsif type == :keyreq
            raise "Value #{p_key} must be set in the environment as #{prefix}#{p_key.to_s.upcase}"
          end
        end
      end

      self.class.new(*new_p_params, **new_k_params).tap do |cursor|
        if (updated_at = value['updated_at'])
          cursor.updated_at = updated_at
        end
      end
    end

    def reset!(&block)
      load_from_db
      @updated_at = Time.at(0)
      block.call unless block.nil?
      save_to_db
      self
    end

    def load_from_db
      if self[:from_env]
        raise "Cursor #{self} was loaded from environment, cannot load from db."
      end

      existing = Polyphemus.instance.db[:cursors].where(name: name).first
      if existing
        @value = existing[:value].to_h
        @updated_at = existing[:updated_at]
        @version = existing[:version]
      end

      self
    end

    def save_to_db
      if self[:from_env]
        raise "Cursor #{self} was loaded from environment, cannot save to db."
      end

      value = Sequel.pg_json_wrap(self.value)
      Polyphemus.instance.db.transaction do
        Polyphemus.instance.db[:cursors].insert(
            name: name,
            updated_at: updated_at,
            value: value,
            version: @version,
        )
      end

      return true
    rescue Sequel::UniqueConstraintViolation => e
      rows_updated = Polyphemus.instance.db[:cursors].where(name: name, version: @version).update(updated_at: updated_at, value: value, version: @version + 1)
      @version += 1
      # If we could not update, it may mean a concurrent modification occurred, so we should load ourselves from the db.
      if rows_updated == 0
        self.load_from_db
        return false
      end

      return true
    end
  end

  class FixedSpanEtlCursor
    attr_reader :name

    def initialize(
      name,
      start_time = ENV['START_TIME'],
      end_time = ENV['END_TIME']
    )
      @name = name
      @value = {
        # 2022-01-18T12:45:57-08:00
        'start_time' => DateTime.parse(start_time),
        'end_time' => DateTime.parse(end_time),
      }
    end

    def reset!(*args, &block)
      raise "#{self.class.name} does not support cursor reset"
    end

    def to_s
      @value.inspect
    end

    def [](k)
      value[k.to_s]
    end

    def []=(k, v)
      value[k.to_s] = v
    end

    def load_from_db
    end

    def save_to_db
    end
  end

  # A group of etl cursors that can sorted and worked upon in 'oldest first' fashion.
  class EtlCursorGroup
    attr_reader :cursors

    def initialize(cursors = [])
      @cursors = cursors
    end

    def with_next(&block)
      sorted_cursors = @cursors.sort_by(&:updated_at)

      sorted_cursors.each do |cursor|
        processed = yield cursor
        if processed
          return true
        end
      end

      false
    end

    def reset_all!
      @cursors.each(&:reset!)
    end
  end
end