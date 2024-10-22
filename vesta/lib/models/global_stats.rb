class Vesta
  class GlobalStats < Sequel::Model(:global_stats)

    def to_hash
      super.to_hash.yield_self do |h|
        h[:recorded_at] = h[:recorded_at].iso8601
        h.delete :id
        h
      end
    end
  end
end
