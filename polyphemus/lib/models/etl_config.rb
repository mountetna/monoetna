class Polyphemus
  class EtlConfig < Sequel::Model
    RUN_ONCE=0
    RUN_NEVER=-1
    STATUS_COMPLETED='completed'
    STATUS_ERROR='error'

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
    plugin :column_encryption do |enc|
      enc.key 0, Polyphemus.instance.config(:etl_secrets_key)
      enc.column :secrets, format: :json
    end

    def etl_job_class
      Polyphemus::Job.from_name(etl)
    end

    def valid_secrets?(secrets)
      (secrets.keys - etl_job_class.secret_keys).empty?
    end

    def run!
      update(ran_at: DateTime.now)

      etl_job_class.new(
        request_params: {},
        request_env: {},
        response: {},
        user: nil
      ).run

      update(
        status: STATUS_COMPLETED,
        output: $stdout.string,
        run_interval: run_interval == Polyphemus::EtlConfig::RUN_ONCE ?
          Polyphemus::EtlConfig::RUN_NEVER : run_interval
      )
    end

    def set_error!(e)
      update(
        status: STATUS_ERROR,
        output: $stdout.string + "\n#{e.message}"
        #run_interval: Polyphemus::EtlConfig::RUN_NEVER
      )
    end

    def as_json
      to_hash.map do |key,value|
        case key
        when :secrets
          next [ key, secrets.map do |skey, secret|
            [ skey, secret.to_s.empty? ? '' : '***' ]
          end.to_h ]
        when :id
          next
        when Time
          next [ key, value.iso8601 ]
        else
          next [ key, value ]
        end
      end.compact.to_h
    end

    def to_revision
      as_json.slice(:config, :updated_at, :comment)
    end

    def validate_config(config)
      schema = JSONSchemer.schema(
        JSON.parse(etl_job_class.as_json[:schema].to_json)
      )
      schema.valid?(JSON.parse(config.to_json))
    end
  end
end
