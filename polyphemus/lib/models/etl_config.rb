class Polyphemus
  class EtlConfig < Sequel::Model
    RUN_ONCE=0
    RUN_NEVER=-1
    RUN_ARCHIVED=-2
    STATUS_COMPLETED='completed'
    STATUS_ERROR='error'

    class << self
      def next_to_run
        self.current.where(
          etl: Polyphemus::Job.list.select(&:should_run?).map(&:job_name)
        ).where {
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
    end

    plugin :timestamps, update_on_create: true, allow_manual_update: true
    plugin :column_encryption do |enc|
      enc.key 0, Polyphemus.instance.config(:etl_secrets_key)
      enc.column :secrets, format: :json
    end

    def self.current
      return self.reverse(:config_id, :version_number).distinct(:config_id).from_self.exclude(run_interval: Polyphemus::EtlConfig::RUN_ARCHIVED)
    end

    def etl_job_class
      Polyphemus::Job.from_name(etl)
    end

    def validate_secrets(secrets)
      etl_job_class.validate_secrets(secrets)
    end

    def validate_params(params)
      etl_job_class.validate_params(params)
    end

    def self.next_id
      Polyphemus.instance.db.get { nextval("etl_configs_ids") }
    end

    def run!
      update(ran_at: DateTime.now)

      task_token = Etna::Clients::Janus.new(
        token: Polyphemus.instance.config(:polyphemus)[:token],
        ignore_ssl: Polyphemus.instance.config(:ignore_ssl),
        **Polyphemus.instance.config(:janus)
      ).generate_token(
        'task',
        signed_nonce: nil,
        project_name: project_name
      )

      job = etl_job_class.new(
        request_params: secrets.merge(params).symbolize_keys.merge(
          project_name: project_name,
          config: config
        ),
        token: task_token
      )

      raise 'Current config is invalid against current schema' unless valid_config?(config)

      job.validate

      raise "Errors in job request: #{job.errors}" unless job.valid?

      results = job.run

      state = {
        status: STATUS_COMPLETED,
        output: $stdout.string
      }

      if run_once?
        state.update(
          run_interval: Polyphemus::EtlConfig::RUN_NEVER,
          params: {}
        )
      end

      update(state)
    end

    def run_once?
      run_interval == Polyphemus::EtlConfig::RUN_ONCE
    end

    def set_error!(e)
      update(
        status: STATUS_ERROR,
        output: $stdout.string + "\n#{e.message}",
        run_interval: Polyphemus::EtlConfig::RUN_NEVER
      )
    end

    def as_json(options = nil)
      to_hash.map do |key,value|
        case key
        when :secrets
          next [ key, secrets.map do |skey, secret|
            [ skey, secret.to_s.empty? ? '' : '***' ]
          end.to_h ]
        when :output, :id
          next
        else
          case value
          when Time
            next [ key, value.iso8601 ]
          end

          next [ key, value ]
        end
      end.compact.to_h
    end

    def to_revision
      as_json.slice(:config, :created_at, :version_number, :comment)
    end

    def valid_config?(config)
      schema = JSONSchemer.schema(
        JSON.parse(etl_job_class.as_json[:schema].to_json)
      )
      schema.valid?(JSON.parse(config.to_json))
    end

    def validate_config(config)
      schema = JSONSchemer.schema(
        JSON.parse(etl_job_class.as_json[:schema].to_json)
      )
      schema.validate(JSON.parse(config.to_json)).to_a
    end
  end
end
