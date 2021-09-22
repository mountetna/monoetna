class Polyphemus
  # This data could eventually get stored in a database, where we track
  #   job status, submission time, user, etc.
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

      def job_name
        self.name.match(/::(.*)Job/)[1].underscore
      end
    end

    attr_reader :request_params, :request_env, :response, :errors, :user
    def initialize(request_params:, request_env:, response:, user:)
      @request_env = request_env
      @request_params = request_params
      @user = user
      @errors = []
      @response = response
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
