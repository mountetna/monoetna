require 'action_controller'

class Polyphemus
  class Controller < Etna::Controller
    def initialize(request, action=nil)
      super
    end
  end

  class StreamingController < ActionController::Base
    include ActionController::Live

    attr_reader :response

    # Need to include these two parameters because our Rack
    #   middleware passes them along when standing up the controller
    #   instance.
    def initialize(request, action=nil)
      super()
      # Copy these over from Etna::Controller because we'll need them.
      @request = request
      @action = action
      # @response = Rack::Response.new
      @params = @request.env['rack.request.params']
      @errors = []
      # @server = @request.env['etna.server']
      # @logger = @request.env['etna.logger']
      @user = @request.env['etna.user']
      # @request_id = @request.env['etna.request_id']
      # @hmac = @request.env['etna.hmac']
    end

    def require_params(*params)
      missing_params = params.reject{|p| @params.key?(p) }
      raise Etna::BadRequest, "Missing param #{missing_params.join(', ')}" unless missing_params.empty?
    end
    alias_method :require_param, :require_params

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
  end
end
