require 'date'
require 'yaml'

class Polyphemus
  # An EtlCursor represents a named "Cursor", generally encoding things such as an updated_at, limit, and offset,
  # which can correctly capture the entire state of query for batched etl jobs.
  class EtlCursor
    attr_reader :value, :name
    attr_accessor :updated_at
    extend ::Etna::Injection::FromHash

    def initialize(name, updated_at = Time.at(0), value = {}, version = 0)
      @name = name
      @value = value
      @updated_at = updated_at
      @version = version
    end

    def load_batch_params(updated_at: nil, batch_end_at: nil)
      return if updated_at.nil? && batch_end_at.nil?
      updated_at = Time.at(0) if updated_at.nil?

      raise "updated_at was not a Time value!" unless updated_at.is_a?(Time)
      raise "batch_end_at was not a Time value!" unless batch_end_at.is_a?(Time)
      @updated_at = updated_at
      self[:batch_end_at] = batch_end_at
    end

    def from_env?
      !self[:batch_end_at].nil?
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

    def reset!(&block)
      load_from_db
      @updated_at = Time.at(0)
      block.call unless block.nil?
      save_to_db
      self
    end

    def load_from_db
      if from_env?
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
      if from_env?
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