# This class handles the http request and routing.
module Etna
  class Server
    class << self
      def route(method, path, options={}, &block)
        @routes ||= []

        @routes << Etna::Route.new(
          method,
          path,
          (@default_options || {}).merge(options),
          &block
        )
      end

      def find_route(request)
        @routes.find do |route|
          route.matches? request
        end
      end

      def with(options={}, &block)
        @default_options = options
        instance_eval(&block)
        @default_options = nil
      end

      def get(path, options={}, &block)
        route('GET', path, options, &block)
      end

      def post(path, options={}, &block)
        route('POST', path, options, &block)
      end

      def put(path, options={}, &block)
        route('PUT', path, options, &block)
      end

      def delete(path, options={}, &block)
        route('DELETE', path, options, &block)
      end

      def route_path(request,name,params={})
        route = routes.find do |route|
          route.name.to_s == name.to_s
        end
        return route ? route.path(params) : nil
      end

      attr_reader :routes
    end

    def call(env)
      request = Rack::Request.new(env)

      request.env['etna.server'] = self

      application.logger.log_request(request)

      route = self.class.find_route(request)

      if route
        @params = request.env['rack.request.params']
        return route.call(self, request)
      end

      [404, {}, ["There is no such path '#{request.path}'"]]
    rescue => e
      application.logger.log_error(e)
      raise
    end

    def initialize
      # Setup logging.
      application.setup_logger

      # This needs to be required before yabeda invocation, but cannot belong at the top of any module since clients
      # do not install yabeda.
      require 'yabeda'
      application.setup_yabeda
    end

    private

    # The base application class is a singleton independent of this rack server,
    # holding e.g. configuration.
    def application
      @application ||= Etna::Application.instance
    end
  end
end
