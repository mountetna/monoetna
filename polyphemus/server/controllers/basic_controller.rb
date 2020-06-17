class BasicController

  def initialize(request, action)

    @request = request
    @params = request.POST()
    @action = action

    @data = nil
  end

  def admin_user?(uri, data)

    begin

      response = JSON.parse(make_request(Secrets::JANUS_ADDR+uri, data))
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
      https_conn.use_ssl = true
      https_conn.verify_mode = OpenSSL::SSL::VERIFY_PEER

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

      raise_err(:SERVER_ERR, 1, m)
    end
  end

  # Proxy the exception.
  def raise_err(type, id, method)

    raise BasicError.new(type, id, method)
  end
end