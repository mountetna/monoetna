class Polyphemus
  # This data could eventually get stored in a database, where we track
  #   job status, submission time, user, etc.
  class Job
    attr_reader :request_params, :request_env, :response, :job_params, :errors, :user
    def initialize(request_params:, request_env:, response:, user:)
      @request_env = request_env
      @request_params = request_params
      @job_params = request_params[:job_params]
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

    def require_job_params(*params)
      missing_params = params.reject{|p| @job_params.key?(p) }
      raise Etna::BadRequest, "job_params missing required param(s): #{missing_params.join(', ')}" unless missing_params.empty?
    end
    alias_method :require_job_param, :require_job_params
  end

  class JobType < String
    REDCAP = JobType.new("redcap")
  end
end
