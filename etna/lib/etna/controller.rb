require 'erb'
require 'net/smtp'

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
    end

    def application
      Etna::Application.instance
    end

    def log(line)
      @logger.warn(request_msg(line))
    end

    def event_log(params)
      begin
        Etna::Application.instance.event_log({
          project_name: @params[:project_name],
          user: @user
        }.compact.merge(params))
      rescue Exception => e
        log("event_log failed with #{e.backtrace} #{e.message}")
      end
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
      return instance_eval(&block) if block_given?
      return send(@action) if @action

      [501, {}, ['This controller is not implemented.']]
    rescue Exception => e
      error = e
    ensure
      log_request if !@request.env['etna.dont_log'] || error
      return handle_error(error) if error
    end

    def try_stream(content_type, &block)
      if @request.env['rack.hijack?']
        @request.env['rack.hijack'].call
        stream = @request.env['rack.hijack_io']

        headers = [
          "HTTP/1.1 200 OK",
          "Content-Type: #{content_type}"
        ]
        stream.write(headers.map { |header| header + "\r\n" }.join)
        stream.write("\r\n")
        stream.flush

        Thread.new do
          block.call(stream)
        ensure
          stream.close
        end

        # IO is now streaming and will be processed by above thread.
        @response.close
      else
        @response['Content-Type'] = content_type
        block.call(@response)
        @response.finish
      end
    end

    def send_email(to_name, to_email, subject, content)
      message = <<MESSAGE_END
From: Data Library <noreply@janus>
To: #{to_name} <#{to_email}>
Subject: #{subject}

#{content}
MESSAGE_END

      unless application.test?
        Net::SMTP.start('smtp.ucsf.edu') do |smtp|
          smtp.send_message message, 'noreply@janus', to_email
        end
      end
    rescue => e
      @logger.log_error(e)
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
      @request.scheme + '://' + application.host + path
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
      [:janus, :magma, :timur, :metis, :vulcan, :polyphemus, :gnomon].map do |host|
        [ :"#{host}_host", application.config(host)&.dig(:host) ]
      end.to_h.compact
    end

    private

    def redact_keys
      @request.env['etna.redact_keys']
    end

    def add_redact_keys(new_redact_keys=[])
      @request.env['etna.redact_keys'] = (@request.env['etna.redact_keys'] || []).concat(new_redact_keys)
    end

    def log_request
      censor = Etna::Censor.new(redact_keys)

      redacted_params = @params.map do |key,value|
        [ key, censor.redact(key, value) ]
      end.to_h

      log("User #{@user ? @user.email : :unknown} calling #{controller_name}##{@action} with params #{redacted_params}")
    end

    def controller_name
      self.class.name.sub("Kernel::", "").sub("Controller", "").downcase
    end

    def success(msg, content_type='text/plain', disposition='not given')
      @response['Content-Type'] = content_type
      @response['Content-Disposition'] = disposition unless disposition=='not given'
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
