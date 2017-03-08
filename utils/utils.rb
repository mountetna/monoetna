require './conf.rb'

module Utils

  def self.make_request(url, data)

    begin

      uri = URI.parse(url)
      https_conn = Net::HTTP.new(uri.host, uri.port)
      #https_conn.use_ssl = true
      #https_conn.verify_mode = OpenSSL::SSL::VERIFY_PEER

      request = Net::HTTP::Post.new(uri.path)
      request.set_form_data(data)

      response = https_conn.request(request)
      return response.body
    rescue Timeout::Error, 
             Errno::EINVAL, 
             Errno::ECONNRESET, 
             EOFError, 
             Net::HTTPBadResponse, 
             Net::HTTPHeaderSyntaxError, 
             Net::ProtocolError => error

      puts 'error: '+url
      return nil
    end
  end

  def self.verify_login(response)

    m = __method__.to_s
    if response.key?('success')

      if response['success']

        if response.key?('user_info')

          if response['user_info'].key?('token')

            puts m+': Saving Token'
            return response['user_info']['token']
          else

            puts m+': Malformed Response'
            return nil
          end
        else

          puts m+': Malformed Response'
          return nil
        end
      end
    else

      puts m+': Malformed Response'
      return nil
    end
    puts Conf::HR
  end
end