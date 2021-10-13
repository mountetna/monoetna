require 'date'

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

    def reset!(&block)
      load_from_db
      @updated_at = Time.at(0)
      block.call unless block.nil?
      save_to_db
      self
    end

    def load_from_db
      existing = Polyphemus.instance.db[:cursors].where(name: name).first
      if existing
        @value = existing[:value].to_h
        @updated_at = existing[:updated_at]
        @version = existing[:version]
      end

      self
    end

    def save_to_db
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