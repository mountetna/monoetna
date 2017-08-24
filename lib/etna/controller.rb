module Etna
  class Controller
    def initialize(request, action=nil)
      @request = request
      @action = action
      @response = Rack::Response.new
      @params = @request.env['rack.request.params']
      @errors = []
      @logger = @request.env['rack.logger']
    end

    def log(line)
      @logger.write "#{Time.now.strftime("%d/%b/%Y %H:%M:%S")} #{line}\n"
    end

    def response
      return send(@action) if @action 

      [501, {}, ['This controller is not implemented.']]
    rescue Etna::BadRequest, Etna::ServerError => e
      return failure(e.status, e.message)
    end

    VIEW_PATH=File.expand_path("../views", __dir__)

    def view(name)
      txt = File.read("#{self.class.VIEW_PATH}/#{name}.html")
      @response['Content-Type'] = 'text/html'
      @response.write(txt)
      @response.finish
    end

    def erb_view(name)
      txt = File.read("#{self.class.VIEW_PATH}/#{name}.html.erb")
      @response['Content-Type'] = 'text/html'
      @response.write(ERB.new(txt).result)
      @response.finish
    end

    private

    def success(content_type, msg)
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
