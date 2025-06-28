class Polyphemus
  class RuntimeConfig < Sequel::Model(:runtime_configs)
    plugin :timestamps, update_on_create: true
    include WithLogger

    # TODO: maybe move the scheduling logic to CronWorkflows
    def self.eligible_runtime_configs
      # Get all runtime configs where interval is not null/none - runtime
      # configs initiated as empty when configs are created
      self.exclude(disabled: true)
        .exclude(run_interval: nil)
        .exclude(run_interval: 0)
        .all.select(&:should_run?)
    end

    def workflow_config
      Polyphemus::Config.current.where(config_id: config_id).first
    end

    def self.for_config(config_id)
      Polyphemus::RuntimeConfig.where(
        config_id: config_id
      ).first
    end

    def workflow_config
      @workflow_config ||= Polyphemus::Config.current.where(
        config_id: config_id,
      ).first
    end

    def last_run
      Polyphemus::Run.where(
        config_id: config_id,
      ).order(:created_at).last
    end
    
    def disabled?
      disabled || run_interval.nil? || run_interval <= 0
    end

    def should_run?
      return false if disabled?

      run = last_run
      
      return true if !run

      if run.is_finished? && run.is_succeeded?
        time_since_last_run = Time.now - DateTime.parse(run.finished_at).to_time

        return time_since_last_run >= run_interval
      end
    end
  end
end
