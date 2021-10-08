require 'erb'

module Etna
  class Controller
    def initialize(request, action = nil)
      @request = request
      @action = action
      @response = Rack::Response.new
      @params = @request.env['rack.request.params']
      @errors = []
      @server = @request.env['etna.server']
      @logger = @request.env['etna.logger']
      @user = @request.env['etna.user']
      @request_id = @request.env['etna.request_id']
      @hmac = @request.env['etna.hmac']

      @route_log_redact_keys = @request.env['etna.route_log_redact_keys']
      @censor = Etna::Censor.new(redact_keys)
    end

    def log(line)
      @logger.warn(request_msg(line))
    end

    def handle_error(e)
      case e
      when Etna::Error
        Rollbar.error(e)
        @logger.error(request_msg("Exiting with #{e.status}, #{e.message}"))
        return failure(e.status, error: e.message)
      else
        Rollbar.error(e)
        @logger.error(request_msg('Caught unspecified error'))
        @logger.error(request_msg(e.message))
        e.backtrace.each do |trace|
          @logger.error(request_msg(trace))
        end
        return failure(500, error: 'Server error.')
      end
    end

    def response(&block)
      log_request

      return instance_eval(&block) if block_given?
      return send(@action) if @action

      [501, {}, ['This controller is not implemented.']]
    rescue Exception => e
      handle_error(e)
    end

    def require_params(*params)
      missing_params = params.reject{|p| @params.key?(p) }
      raise Etna::BadRequest, "Missing param #{missing_params.join(', ')}" unless missing_params.empty?
    end
    alias_method :require_param, :require_params

    def route_path(name, params={})
      @server.class.route_path(@request, name, params)
    end

    def route_url(name, params={})
      path = route_path(name,params)
      return nil unless path
      @request.scheme + '://' + @request.host + path
    end

    # methods for returning a view
    VIEW_PATH = :VIEW_PATH

    def view(name)
      txt = File.read("#{self.class::VIEW_PATH}/#{name}.html")
      @response['Content-Type'] = 'text/html'
      @response.write(txt)
      @response.finish
    end

    def erb_partial(name)
      txt = File.read("#{self.class::VIEW_PATH}/#{name}.html.erb")
      ERB.new(txt).result(binding)
    end

    def erb_view(name)
      @response['Content-Type'] = 'text/html'
      @response.write(erb_partial(name))
      @response.finish
    end

    def config_hosts
      [:janus, :magma, :timur, :metis, :vulcan, :polyphemus].map do |host|
        [ :"#{host}_host", @server.send(:application).config(host)&.dig(:host) ]
      end.to_h.compact
    end

    private

    def redact_keys
      # Subclasses may want to override this, if they want to
      #   redact additional keys
      @route_log_redact_keys
    end

    def log_request
      redacted_params = @params.map do |key,value|
        [ key, @censor.redact(key, value) ]
      end.to_h

      log("User #{@user ? @user.email : :unknown} calling #{controller_name}##{@action} with params #{redacted_params}")
    end

    def controller_name
      self.class.name.sub("Kernel::", "").sub("Controller", "").downcase
    end

    def success(msg, content_type='text/plain')
      @response['Content-Type'] = content_type
      @response.write(msg)
      @response.finish
    end

    def success_json(params)
      success(params.to_json, 'application/json')
    end

    def failure(status, msg)
      @response['Content-Type'] = 'application/json'
      @response.status = status
      @response.write(msg.to_json)
      @response.finish
    end

    def success?
      @errors.empty?
    end

    def error(msg)
      if msg.is_a?(Array)
        @errors.concat(msg)
      else
        @errors.push(msg)
      end
    end

    def request_msg(msg)
      "#{@request_id} #{msg}"
    end
  end
end
