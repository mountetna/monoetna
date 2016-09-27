# metis.rb

# This class handles the http request and routing
class Metis

  def initialize()

    @routes = {}
    @request = {}
    @redis_service = RedisService.new()
  end

  def call(env)

    # Parse the request
    @request = Rack::Request.new(env)
    route = @routes[[@request.request_method, @request.path]]

    if route

      call_action_for(route)
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
    controller_class.new(@redis_service, @request, action).run()
  end
end