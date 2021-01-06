require_relative 'controller'
require_relative '../models/redcap_job'
require_relative '../../etls/redcap/redcap_etl_script_runner'


class JobController < Polyphemus::Controller
  def submit
    require_params(:project_name, :job_type, :job_params)

    case @params[:job_type]
    when Polyphemus::JobType::REDCAP
      job = Polyphemus::RedcapJob.new(@params, @request.env, @response, @user)
    end

    raise Etna::BadRequest, "Unsupported job type: #{@params[:job_type]}." unless job

    job.validate

    raise Etna::BadRequest, "Errors in job request: #{job.errors}" unless job.valid?

    results = job.run

    success({results: results}.to_json, 'application/json')
  rescue => e
    Polyphemus.instance.logger.log_error(e)
    return failure(422, e.message)
  end
end
