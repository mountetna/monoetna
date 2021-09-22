require_relative 'controller'
require_relative '../../etls/redcap/redcap_etl_script_runner'

class JobController < Polyphemus::Controller
  def submit
    require_params(:project_name, :job_type, :job_params)

    job_type = "#{@params[:job_type]}".to_s.split('_').map(&:capitalize).join
    class_name = "Polyphemus::#{job_type}Job"
    job_class = Kernel.const_defined?(class_name) ?  Kernel.const_get(class_name) : nil

    raise "Unsupported job type: #{@params[:job_type]}." unless job_class

    job = job_class.new(
      request_params: @params,
      request_env: @request.env,
      response: @response,
      user: @user)

    job.validate

    raise Etna::BadRequest, "Errors in job request: #{job.errors}" unless job.valid?

    results = job.run

    success({results: results}.to_json, 'application/json')
  rescue => e
    Polyphemus.instance.logger.log_error(e)
    return failure(422, e.message)
  end
end
