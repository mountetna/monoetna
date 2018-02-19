# metis.rb
require 'sequel'
require 'extlib'

# This class handles the http request and routing
class Metis
  include Etna::Application

  attr_reader :db

  def setup_db
    @db = Sequel.connect(config(:db))
    @db.extension :connection_validator
    @db.pool.connection_validation_timeout = -1
  end

  def load_models
    setup_db

    require_relative 'models'
  end

  def project_path(project)
    config(:project_paths)[project.to_sym]
  end

  # Routes are added in the './routes.rb' file
  def add_route(method, path, handler)
    @routes[[method, path]] = handler
  end

  private 
  def call_action_for(route)

    controller, action = route.split('#')
    controller_class = Kernel.const_get(controller)
    controller_class.new(@request, action).run()
  end

  def send_err(err)

    ip = @request.env['HTTP_X_FORWARDED_FOR'].to_s
    ref_id = SecureRandom.hex(4).to_s
    response = { :success=> false, :ref=> ref_id }
    m = err.method.to_s

    case err.type
    when :SERVER_ERR

      code = Conf::ERRORS[err.id].to_s
      @app_logger.error(ref_id+' - '+code+', '+m+', '+ip)
      response[:error] = 'Server error.'
    when :BAD_REQ

      code = Conf::WARNS[err.id].to_s
      @app_logger.warn(ref_id+' - '+code+', '+m+', '+ip)
      response[:error] = 'Bad request.'
    when :BAD_LOG

      code = Conf::WARNS[err.id].to_s
      @app_logger.warn(ref_id+' - '+code+', '+m+', '+ip)
      response[:error] = 'Invalid login.'
    else

      @app_logger.error(ref_id+' - UNKNOWN, '+m+', '+ip)
      response[:error] = 'Unknown error.'
    end

    return response
  end
end
