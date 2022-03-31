require 'json'
require_relative './vulcan_controller'

class SessionsController < Vulcan::Controller
  def session
    begin
      @the_session ||= Session.from_json(@params.slice(:key, :project_name, :workflow_name, :inputs))
    rescue => e
      raise Etna::BadRequest.new(e.message)
    end
  end

  def orchestration
    session.orchestration
  end

  def scheduler
    orchestration.scheduler
  end

  def workflow
    orchestration.workflow
  end

  def submit
    session.run!(storage: storage, token: create_task_token)
    success_json(session.status(storage: storage))
  rescue => e
    Vulcan.instance.logger.log_error(e)
    raise Etna::BadRequest.new(e.message)
  end

  def status
    success_json(session.status(storage: storage))
  rescue => e
    Vulcan.instance.logger.log_error(e)
    raise Etna::BadRequest.new(e.message)
  end

  def create_task_token
    # For development, we can inject a production token
    #   via config.yml, to talk directly to production services.
    #   This should reduce duplication of data.
    # Note: You'll need to configure the Magma host
    #   in config.yml :development to also point to
    #   production Magma.
    if :development == Vulcan.instance.environment
      archimedes_token = Vulcan.instance.config(:archimedes_token)
      return archimedes_token unless archimedes_token.nil?
    end

    janus_client = Etna::Clients::Janus.new(token: @user.token, host: Vulcan.instance.config(:janus)[:host])
    janus_client.generate_token('task', project_name: @params[:project_name], read_only: true)
  end
end

