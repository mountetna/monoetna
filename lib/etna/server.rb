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
      @request = Rack::Request.new(env)

      @request.env['etna.server'] = self

      @logger.request = @request

      route = self.class.routes.find do |route|
        route.matches? @request
      end

      if route
        @params = @request.env['rack.request.params']
        return route.call(self, @request)
      end

      [404, {}, ["There is no such path '#{@request.path}'"]]
    end

    def initialize(config)

      # Setup the application, since we are booting through the rack server.
      application.configure(config)

      # Setup logging.
      setup_logger
    end

    private

    # The base application class is a singleton independent of this rack server,
    # holding e.g. configuration.
    def application
      @application ||= Etna::Application.instance
    end

    def setup_logger
      @logger = Etna::Logger.new(
        # The name of the log_file, required.
        application.config(:log_file),
        # Number of old copies of the log to keep.
        application.config(:log_copies) || 5,
        # How large the log can get before overturning.
        application.config(:log_size) || 1048576
      )
      log_level = (application.config(:log_level) || 'warn').upcase.to_sym

      @logger.level = Logger.const_defined?(log_level) ? Logger.const_get(log_level) : Logger::WARN

      # the application logger is available to anyone
      application.set_logger(@logger)
    end
  end
end
