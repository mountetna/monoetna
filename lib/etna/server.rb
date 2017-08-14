# This class handles the http request and routing
module Etna
  class Server
    class << self
      def route path, method=nil, action=nil, &block
        @routes ||= {}

        @routes[path] = Etna::Route.new(
          method,
          action,
          &block
        )
      end

      def get path, action=nil, &block
        route(path, 'GET', action, &block)
      end
      
      def post path, action=nil, &block
        route(path, 'POST', action, &block)
      end

      def put path, action=nil, &block
        route(path, 'PUT', action, &block)
      end

      def delete path, action=nil, &block
        route(path, 'DELETE', action, &block)
      end

      attr_reader :routes
    end

    def call(env)
      @request = Rack::Request.new(env)
      route = @routes[[@request.request_method, @request.path]]

      if self.class.routes.has_key? @request.path
        return self.class.routes[@request.path].call(self, @request)
      end

      [ 404, {}, ["There is no such path '#{@request.path}'"] ]
    end
  end
end
