class BasicController

  def initialize(request, action, logger)

    @request = request
    @params = request.POST()
    @action = action
    @logger = logger
  end

  def admin_user?(uri, data)

    begin

      response = JSON.parse(make_request(Conf::JANUS_ADDR+uri, data))
      if response.key?('administrator') && response['administrator']

        return true
      else

        return false
      end
    rescue

      return false
    end
  end

  def make_request(url, data)

    m = __method__
    begin

      uri = URI.parse(url)
      https_conn = Net::HTTP.new(uri.host, uri.port)
      #https_conn.use_ssl = true
      #https_conn.verify_mode = OpenSSL::SSL::VERIFY_PEER

      request = Net::HTTP::Post.new(uri.path)
      request.set_form_data(data)

      response = https_conn.request(request)
      response_code = response.code.to_i()

      return (response_code == 200) ? response.body : send_err(:SERVER_ERR,0,m)
    rescue

      send_err(:SERVER_ERR, 1, m)
    end
  end

  def send_err(type, id, method)

    ip = @request.env['HTTP_X_FORWARDED_FOR'].to_s
    ref_id = SecureRandom.hex(4).to_s
    response = { :success=> false, :ref=> ref_id }

    case type
    when :SERVER_ERR

      code = Conf::ERRORS[id].to_s
      @logger.error(ref_id+' - '+code+', '+method.to_s+', '+ip)
      response[:error] = 'Server error.'
    when :BAD_REQ

      code = Conf::WARNS[id].to_s
      @logger.warn(ref_id+' - '+code+', '+method.to_s+', '+ip)
      response[:error] = 'Bad request.'
    when :BAD_LOG

      code = Conf::WARNS[id].to_s
      @logger.warn(ref_id+' - '+code+', '+method.to_s+', '+ip)
      response[:error] = 'Invalid login.'
    else

      @logger.error(ref_id+' - UNKNOWN, '+method.to_s+', '+ip)
      response[:error] = 'Unknown error.'
    end

    return response.to_json()
  end
end