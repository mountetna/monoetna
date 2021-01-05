require 'action_controller'

class Polyphemus
  class Controller < Etna::Controller
    def initialize(request, action=nil)
      super
    end
  end

  class StreamingController < Etna::Controller
    # include ActionController::Live


    # Need to include these two parameters because our Rack
    #   middleware passes them along when standing up the controller
    #   instance.
    # def initialize(request, action=nil)
    #   super()

    #   # Duplicate from Etna::Controller
    #   @request = request
    #   @action = action
    #   @params = @request.env['rack.request.params']
    #   @errors = []
    #   @server = @request.env['etna.server']
    #   @logger = @request.env['etna.logger']
    #   @user = @request.env['etna.user']
    #   @request_id = @request.env['etna.request_id']
    #   @hmac = @request.env['etna.hmac']

    #   # Override this, to make response an ActionController::Live streaming response.
    #   # @response = ActionController::Live.make_response!(request)
    # end

    # def require_params(*params)
    #   missing_params = params.reject{|p| @params.key?(p) }
    #   raise Etna::BadRequest, "Missing param #{missing_params.join(', ')}" unless missing_params.empty?
    # end
    # alias_method :require_param, :require_params

    # def response(&block)
    #   return instance_eval(&block) if block_given?

    #   return send(@action) if @action


    #   [501, {}, ['This controller is not implemented.']]
    # rescue Etna::Error => e
    #   Rollbar.error(e)
    #   @logger.error(request_msg("Exiting with #{e.status}, #{e.message}"))
    #   return failure(e.status, error: e.message)
    # rescue Exception => e
    #   Rollbar.error(e)
    #   @logger.error(request_msg('Caught unspecified error'))
    #   @logger.error(request_msg(e.message))
    #   e.backtrace.each do |trace|
    #     @logger.error(request_msg(trace))
    #   end
    #   return failure(500, error: 'Server error.')
    # end

    # private

    # def success(msg, content_type='text/plain')
    #   @response['Content-Type'] = content_type
    #   @response.write(msg)
    #   @response.finish
    # end

    # def success_json(params)
    #   success(params.to_json, 'application/json')
    # end

    # def failure(status, msg)
    #   @response['Content-Type'] = 'application/json'
    #   @response.status = status
    #   @response.write(msg.to_json)
    #   @response.finish
    # end

    # def success?
    #   @errors.empty?
    # end

    # def error(msg)
    #   if msg.is_a?(Array)
    #     @errors.concat(msg)
    #   else
    #     @errors.push(msg)
    #   end
    # end

    # def request_msg(msg)
    #   "#{@request_id} #{msg}"
    # end
  end
end
