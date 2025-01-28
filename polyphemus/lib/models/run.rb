class Polyphemus
  class Run < Sequel::Model(:runs)
    plugin :timestamps, update_on_create: true

    def started_at
      return nil unless self.orchestrator_metadata
      self.orchestrator_metadata['startedAt']
    end

    def finished_at
      return nil unless self.orchestrator_metadata
      self.orchestrator_metadata['finishedAt']
    end

    def status
      return nil unless self.orchestrator_metadata
      self.orchestrator_metadata['phase']
    end

    def is_finished?
      self.finished_at.present?
    end

    def is_succeeded?
      self.status == 'Succeeded'
    end

    def self.last_for_config(config_id)
      Polyphemus::Run.where(
        config_id: config_id
      ).reverse(:created_at).first
    end

    def summary_hash
      {
        run_id: run_id,
        version_number: version_number,
        config_id: config_id,
        status: status,
        created_at: created_at,
        finished_at: finished_at
      }
    end
  end
end
