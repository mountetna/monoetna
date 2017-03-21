class BasicController

  def initialize(redis_service, request, action)

    @redis_service = redis_service
    @params = request.POST()
    @action = action
    @user = nil
    @data = nil
  end

  def set_user()

    # Check for the user token.
    if !@params.key?('token') then raise_err(:BAD_REQ, 0, __method__) end

    # Get the user and their permissions.
    data = { :token=> @params['token'], :app_key=> Secrets::APP_KEY }
    response = JSON.parse(make_request(Conf::JANUS_ADDR+'/check', data))

    # Check that the user query is valid. 
    if !response.key?('success') then raise_err(:BAD_REQ, 0, __method__) end
    if !response['success'] then raise_err(:BAD_REQ, 2, __method__) end

    @user = UserModel.new(response)
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

      if response_code == 200

        return response.body
      else

        raise_err(:SERVER_ERR,0,m)
      end
    rescue

      raise_err(:SERVER_ERR, 6, m)
    end
  end

  # Proxy the exception.
  def raise_err(type, id, method)

    raise BasicError.new(type, id, method)
  end
end