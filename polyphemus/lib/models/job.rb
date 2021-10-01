class Polyphemus
  # This data could eventually get stored in a database, where we track
  #   job status, submission time, etc.
  class Job
    class << self
      def inherited(subclass)
        @list ||= []
        @list << subclass
      end
      attr_reader :list

      def from_name(job_name)
        list.find do |job|
          job.job_name == job_name
        end
      end

      def secret_keys
        as_json[:secrets]
      end

      def validate_secrets(secrets)
        unless (secrets.keys - secret_keys).empty?
          return "Secrets for #{job_name} jobs must be one of: #{secret_keys.join(', ')}"
        end
      end

      def job_params
        as_json[:params]
      end

      def validate_params(params)
        errors = []
        params.each do |param, value|
          unless job_params.has_key?(param)
            errors.push("no such param #{param}")
            next
          end
          value_opts = job_params[param]
          case value_opts
          when Array
            unless value_opts.include?(value)
              errors.push("#{param} must be in: #{value_opts.join(', ')}")
            end
          when 'string'
            unless value.is_a?(String)
              errors.push("#{param} must be a string")
            end
          when 'boolean'
            unless [true,false].include?(value)
              errors.push("#{param} must be a boolean")
            end
          end
        end
        errors
      end

      def job_name
        self.name.match(/::(.*)Job/)[1].underscore
      end
    end

    attr_reader :request_params, :errors, :token
    def initialize(request_params:, token:)
      @request_params = request_params
      @token = token
      @errors = []
    end

    def run
      raise "This should be implemented in a subclass."
    end

    def valid?
      errors.length == 0
    end

    def validate
      raise "This should be implemented in a subclass."
    end

    private

    def require_params(*params)
      missing_params = params.reject{|p| @request_params.key?(p) }
      raise JobError, "request_params missing required param(s): #{missing_params.join(', ')}" unless missing_params.empty?
    end
    alias_method :require_param, :require_params
  end

  class JobType < String
    REDCAP = JobType.new("redcap")
  end

  class JobError < StandardError
  end
end

require_relative './redcap_job'
