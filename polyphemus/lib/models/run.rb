class Polyphemus
  class Run < Sequel::Model(:runs)

    def started_at
      self.orchestrator_metadata['startedAt']
    end

    def finished_at
      self.orchestrator_metadata['finishedAt']
    end

    def status
      self.orchestrator_metadata['phase']
    end

    def is_finished?
      self.finished_at.present?
    end

    def is_succeeded?
      self.status == 'Succeeded'
    end

  end
end
  