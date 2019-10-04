module Etna
  class DescribeRoutes
    def initialize(app)
      @app = app
    end

    def call(env)
      request = Rack::Request.new(env)

      return @app.call(env) unless request.request_method == 'OPTIONS' && request.path == '/'

      return [ 200, { 'Content-Type' => 'application/json' }, [ route_json ] ]
    end

    private

    def route_json
      server.routes.map do |route|
        route.to_hash
      end.to_json
    end

    def application
      @application ||= Etna::Application.instance
    end

    def server
      @server ||= application.class.const_get(:Server)
    end
  end
end
