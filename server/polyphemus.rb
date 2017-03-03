# This class handles the http request and routing
class Polyphemus

  def initialize()

    @routes = {}
    @request = {}

    # Log file details
    path = ::File.dirname(::File.expand_path(__FILE__))
    log_file = ::File.join(path,'..','log','app.log')
    @app_logger = ::Logger.new(log_file, 5, 1048576)
    @app_logger.level = Logger::WARN
  end

  def call(env)

    # Parse the request
    @request = Rack::Request.new(env)
    route = @routes[[@request.request_method, @request.path]]

    if route

      Rack::Response.new(call_action_for(route))
    else

      Rack::Response.new('File not found.', 404)
    end
  end

  # Routes are added in the './routes.rb' file
  def add_route(method, path, handler)

    @routes[[method, path]] = handler
  end

  private 
  def call_action_for(route)

    controller, action = route.split('#')
    controller_class = Kernel.const_get(controller)
    controller_class.new(@request, action, @app_logger).run()
  end
end