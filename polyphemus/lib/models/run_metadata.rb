class Polyphemus
  class RunMetadata < Sequel::Model(:run_metadata)
    def self.new_entries
      self.where(Sequel.lit('scheduled_at <= ?', Time.now)).all
    end
    def update_scheduled_at(time)
      self.update(scheduled_at: time)
    end
  end
end
