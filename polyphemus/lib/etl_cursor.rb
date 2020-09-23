class Polyphemus
  # An EtlCursor represents a named "Cursor", generally encoding things such as an updated_at, limit, and offset,
  # which can correctly capture the entire state of query for batched etl jobs.
  class EtlCursor
    attr_reader :value, :updated_at

    def initialize(name, updated_at = Time.now, value = {})
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