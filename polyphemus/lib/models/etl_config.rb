class Polyphemus
  class EtlConfig < Sequel::Model
    RUN_ONCE=0
    RUN_NEVER=-1
    STATUS_COMPLETED='completed'

    def self.next_to_run
      self.exclude(archived: true).where {
        Sequel.|(
          {run_interval: Polyphemus::EtlConfig::RUN_ONCE},
          Sequel.&(
            (run_interval > 0),
            Sequel.|(
              {ran_at: nil},
              ((ran_at + Sequel.lit("? * interval '1 second'", :run_interval)) < DateTime.now)
            )
          )
        )
      }.first
    end
    plugin :timestamps, update_on_create: true

    def etl_job_class
      Polyphemus::Job.from_name(etl)
    end

    def run!
      update(ran_at: DateTime.now)

      output = StringIO.new

      old_stdout = $stdout
      $stdout = output

      etl_job_class.new(
        request_params: {},
        request_env: {},
        response: {},
        user: nil
      ).run

      $stdout = old_stdout
      update(
        status: STATUS_COMPLETED,
        output: output.string,
        run_interval: run_interval == Polyphemus::EtlConfig::RUN_ONCE ?
          Polyphemus::EtlConfig::RUN_NEVER : run_interval
      )
    end

    def as_json
      to_hash.reject { |k,v| k == :id }
    end

    def validate_config(config)
      schema = JSONSchemer.schema(
        JSON.parse(etl_job_class.as_json[:schema].to_json)
      )
      schema.valid?(JSON.parse(config.to_json))
    end
  end
end
