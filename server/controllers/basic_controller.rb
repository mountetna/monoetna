class BasicController

  def initialize(request, action)

    @request = request
    @params = request.POST()
    @action = action

    @user = nil
    @file = nil    # postgres model
    @upload = nil  # postgres model
  end

  def set_user()

    # Check for the user token.
    raise_err(:BAD_REQ, 0, __method__) if !@params.key?('token')

    # Get the user and their permissions.
    data = { :token=> @params['token'], :app_key=> Secrets::APP_KEY }
    response = JSON.parse(make_request(Conf::JANUS_ADDR+'/check', data))

    # Check that the user query is valid.
    raise_err(:BAD_REQ, 0, __method__) if !response.key?('success')
    raise_err(:BAD_REQ, 2, __method__) if !response['success']

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

    begin

      uri = URI.parse(url)
      https_conn = Net::HTTP.new(uri.host, uri.port)
      https_conn.use_ssl = true
      https_conn.verify_mode = OpenSSL::SSL::VERIFY_PEER

      request = Net::HTTP::Post.new(uri.path)
      request.set_form_data(data)

      response = https_conn.request(request)
      response_code = response.code.to_i()

      if response_code == 200

        return response.body
      else

        raise_err(:SERVER_ERR, 0, __method__)
      end
    rescue

      raise_err(:SERVER_ERR, 6, __method__)
    end
  end

  # Proxy the exception.
  def raise_err(type, id, method)

    raise BasicError.new(type, id, method)
  end
end