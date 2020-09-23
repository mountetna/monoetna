require 'date'

class Polyphemus
  # An EtlCursor represents a named "Cursor", generally encoding things such as an updated_at, limit, and offset,
  # which can correctly capture the entire state of query for batched etl jobs.
  class EtlCursor
    attr_reader :value, :name
    attr_accessor :updated_at

    def initialize(name, updated_at = nil, value = {})
      @name = name
      @value = value
      @updated_at = updated_at
    end

    def [](k)
      value[k.to_s]
    end

    def []=(k, v)
      value[k.to_s] = v
    end

    def load_from_db
      existing = Polyphemus.instance.db[:cursors].where(name: name).first
      if existing
        @value = existing[:value]
        @updated_at = existing[:updated_at]
      else
        @value = {}
        @updated_at = nil
      end

      self
    end

    def save_to_db
      value = Sequel.pg_json_wrap(self.value)
      Polyphemus.instance.db[:cursors].insert(
          name: name,
          updated_at: updated_at,
          value: value,
      )
    rescue  Sequel::UniqueConstraintViolation => e
      Polyphemus.instance.db[:cursors].where(name: name).update(updated_at: updated_at, value: value)
    end
  end

  # A group of etl cursors that can sorted and worked upon in 'oldest first' fashion.
  class EtlCursorGroup
    def initialize(cursors = [])
      @cursors = cursors
    end

    def with_next(&block)
      n = @cursors.inject(nil) do |acc, n|
        if acc.nil?
          n
        elsif n.updated_at < acc.updated_at
          n
        else
          acc
        end
      end

      yield n unless n.nil?
    end
  end
end