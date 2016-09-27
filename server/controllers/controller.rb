# controller.rb
# The generic controller.

class Controller

  def initialize(redis_service, request, action)

    @redis_service = redis_service
    @request = request
    @action = action
  end

  def run()  

    if !request_valid?()

      send_bad_request()
    else

      send(@action)
    end
  end

  def request_valid?()

    if !@request.post?()

      # POST params are not present.
      return false
    end

    if !SignService::verify_request_parameters(@request.POST())

      # POST params are not in the correct format.
      return false
    end

    if generate_signature(@request.POST()) != @request.POST()['signature']

      # The packet doesn't have the correct HMAC signature/hash
      return false
    end

    true
  end

  def send_bad_request()

    Rack::Response.new({ success: false, error: 'Bad request.' }.to_json())
  end
end