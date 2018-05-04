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
    end

    def log(line)
      @logger.write "#{Time.now.strftime("%d/%b/%Y %H:%M:%S")} #{line}\n"
    end

    def response(&block)
      return instance_eval(&block) if block_given?

      return send(@action) if @action 

      [501, {}, ['This controller is not implemented.']]
    rescue Etna::Error => e
      return failure(e.status, error: e.message)
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

    def erb_view(name)
      txt = File.read("#{self.class::VIEW_PATH}/#{name}.html.erb")
      @response['Content-Type'] = 'text/html'
      @response.write(ERB.new(txt).result(binding))
      @response.finish
    end

    private

    def success(msg, content_type='text/plain')
      @response['Content-Type'] = content_type
      @response.write(msg)
      @response.finish
    end

    def failure(status, msg)
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
  end
end
