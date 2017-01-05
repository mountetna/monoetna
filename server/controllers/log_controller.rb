# log_controller.rb
# This controller serves as an intermediary to Janus logging services.

class LogController

  def initialize(redis_service, request, action)

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

        return send_server_error()
      end
    rescue Timeout::Error, 
           Errno::EINVAL, 
           Errno::ECONNRESET, 
           EOFError, 
           Net::HTTPBadResponse, 
           Net::HTTPHeaderSyntaxError, 
           Net::ProtocolError => error

      return send_server_error()
    end
  end

   def send_bad_login()

    Rack::Response.new({ success: false, error: 'Invalid login.' }.to_json())
  end

  def send_bad_request()

    Rack::Response.new({ success: false, error: 'Bad request.' }.to_json())
  end

  def send_server_error()

    error_message = 'There was a server error.'
    Rack::Response.new({ success: false, error: error_message }.to_json())
  end
end