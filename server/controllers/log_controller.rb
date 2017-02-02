# log_controller.rb
# This controller serves as an intermediary to Janus logging services.

class LogController

  def initialize(redis_service, request, action, logger)

    @request = request
    @action = action
  end

  def run()  

    send(@action)
  end

  def log_in() 

    # Get the params out of the POST
    params = @request.POST()

    # Check for the correct parameters.
    if params.key?('email') && params.key?('pass')

      url = Conf::JANUS_ADDR
      url = url + '/login'
      data = { 

        :email=> params['email'], 
        :pass=> params['pass'], 
        :app_key=> Conf::APP_KEY 
      }
      response = make_request(url, data)
      return response
    else

      return send_bad_request()
    end
  end

  def log_out()

    # Get the params out of the POST
    params = @request.POST()

    # Check for the correct parameters.
    if params.key?('email') && params.key?('token')

      url = Conf::JANUS_ADDR
      url = url + '/logout'
      data = { 

        :email=> params['email'], 
        :token=> params['token'], 
        :app_key=> Conf::APP_KEY 
      }
      response = make_request(url, data)
      return response
    else

      return send_bad_request()
    end
  end

  # This is use for external requests. There is another in the generic
  # controller that is used for internal requests.
  def check_log()

    # Get the params out of the POST
    params = @request.POST()

    if params.key?('token')

      url = Conf::JANUS_ADDR
      url = url + '/check'
      data = { :token=> params['token'], :app_key=> Conf::APP_KEY }
      response = make_request(url, data)
      return response
    else

      return send_bad_request()
    end
  end

  def make_request(url, data)

    uri = URI.parse(url)
    http = Net::HTTP.new(uri.host, uri.port)
    request = Net::HTTP::Post.new(uri.request_uri)
    request.set_form_data(data)

    begin

      response = http.request(request)
      response_code = response.code.to_i()

      if response_code == 200

        return Rack::Response.new(response.body)
      else

        met = __method__.to_s + ', '+ url +', '+ response_code.to_s
        return send_server_error(6, met)
      end
    rescue Timeout::Error, 
           Errno::EINVAL, 
           Errno::ECONNRESET, 
           EOFError, 
           Net::HTTPBadResponse, 
           Net::HTTPHeaderSyntaxError, 
           Net::ProtocolError => error

      met = __method__.to_s + ', '+ url +', '+ response_code.to_s
      return send_server_error(6, met)
    end
  end

   def send_bad_login()

    Rack::Response.new({ success: false, error: 'Invalid login.' }.to_json())
  end

  def send_bad_request(id, method)

    ref_id = SecureRandom.hex(4)
    code = Conf::WARNS[id].to_s
    @logger.warn(ref_id.to_s+' - '+code+', '+method.to_s)
    response = {

      :success=> false,
      :error=> 'Bad request.',
      :reference_id=> ref_id,
      :error_code=> id
    }
    Rack::Response.new(response.to_json())
  end

  def send_server_error(id, method)

    ref_id = SecureRandom.hex(4)
    code = Conf::ERRORS[id].to_s
    @logger.error(ref_id.to_s+' - '+code+', '+method.to_s)
    response = { 

      :success=> false,
      :error=> 'There was a server error.',
      :reference_id=> ref_id,
      :error_code=> id
    }
    Rack::Response.new(response.to_json())
  end
end