class Polyphemus
  class RuntimeConfig < Sequel::Model(:runtime_configs)
    plugin :timestamps, update_on_create: true
    include WithLogger

    # TODO: maybe move the scheduling logic to CronWorkflows
    def self.eligible_runtime_configs

      eligible_runtime_configs = []

      # Get all runtime configs where interval is not null/none - runtime configs initiated as empty when configs are created
      runtime_configs = self.exclude(disabled: true).exclude(run_interval: nil).all

      runtime_configs.select do |runtime_config|
     
        config = Polyphemus::Config.current.where(
          config_id: runtime_config.config_id,
        ).first

        run = Polyphemus::Run.where(
            config_id: config.config_id,
        ).order(:created_at).last
        
        next unless run&.orchestrator_metadata
        
        if run.is_finished?
          if run.is_succeeded?
            # Check if enough time has elapsed since last run based on interval
            time_since_last_run = Time.now - DateTime.parse(run.finished_at).to_time
            if time_since_last_run >= runtime_config.run_interval
              eligible_runtime_configs << runtime_config
            end
          end
        end
      end
      eligible_runtime_configs
    end

    def self.for_config(config_id)
      Polyphemus::RuntimeConfig.where(
        config_id: config_id
      ).first || nil
    end
  end
end
