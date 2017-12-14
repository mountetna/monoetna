# This class handles the http request and routing.
module Etna
  class Server
    class << self
      def route(method, path, options={}, &block)
        @routes ||= []

        @routes << Etna::Route.new(
          method,
          path,
          options,
          &block
        )
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

      attr_reader :routes
    end

    def call(env)
      @request = Rack::Request.new(env)

      @request.env['etna.server'] = self
      @request.env['etna.logger'] = @logger

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
      @application ||= Etna::Application.find(self.class)
    end

    def setup_logger
      @logger = Logger.new(
        # The name of the log_file, required.
        application.config(:log_file),
        # Number of old copies of the log to keep.
        application.config(:log_copies) || 5,
        # How large the log can get before overturning.
        application.config(:log_size) || 1048576
      )
      @logger.level = Logger::WARN
    end
  end
end
