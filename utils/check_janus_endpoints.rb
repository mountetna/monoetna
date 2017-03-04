  # This will check the Janus URL endpoints.
require 'net/http'
require 'openssl'
require 'json'

require '../server/secrets'

class JanusCheck

  def initialize()

    @good_email = 'jasondcater@gmail.com'
    @good_passwd = '123abc'
    @good_token = nil
    @hr = '-----------------------------------------------------------'
  end

  def log_in() 

    url = 'http://janus-dev.ucsf.edu/login'
    data = {

      :email=> @good_email, 
      :pass=> @good_passwd, 
      :app_key=> Secrets::APP_KEY
    }
    response = JSON.parse(make_request(url, data))
    verify_login(response)
    puts response
    puts ''
  end

  def verify_login(response)

    m = __method__.to_s
    if response.key?('success')

      if response['success']

        if response.key?('user_info')

          if response['user_info'].key?('token')

            puts m+': Saving Token'
            @good_token = response['user_info']['token']
          else

            puts m+': Malformed Response'
          end
        else

          puts m+': Malformed Response'
        end
      end
    else

      puts m+': Malformed Response'
    end
    puts @hr
  end

  def check_log()

    m = __method__.to_s
    if !@good_token

      puts m+': Bad Token'
      return
    end

    url = 'http://janus-dev.ucsf.edu/check'
    data = { :token=> @good_token, :app_key=> Secrets::APP_KEY }
    response = JSON.parse(make_request(url, data))
    verify_check_log(response)
    puts response
    puts ''
  end

  def verify_check_log(response)

    m = __method__.to_s
    if response.key?('success')

      if response['success']

        if response.key?('user_info')

          if response['user_info'].key?('permissions')

            puts m+': Permissions Key Checked'
          else

            puts m+': Malformed Response'
          end
        else

          puts m+': Malformed Response'
        end
      end
    else

      puts m+': Malformed Response'
    end
    puts @hr
  end

  def log_out() 

    m = __method__.to_s
    if !@good_token

      puts m+': Bad Token'
      return
    end

    url = 'http://janus-dev.ucsf.edu/logout'
    data = { :token=> @good_token, :app_key=> Secrets::APP_KEY }
    verify_log_out(make_request(url, data))
  end

  def verify_log_out(response)

    puts __method__.to_s+': You should see \'success\'=>false'
    check_log()
  end

  def check_admin()

    url = 'http://janus-dev.ucsf.edu/check-admin'
    data = {

      :email=> @good_email, 
      :pass=> @good_passwd, 
      :app_key=> Secrets::APP_KEY
    }
    response = JSON.parse(make_request(url, data))
    puts response
    puts ''
  end

  def check_admin_token()

    log_in()
    url = 'http://janus-dev.ucsf.edu/check-admin-token'
    data = { :token=> @good_token, :app_key=> Secrets::APP_KEY }
    response = JSON.parse(make_request(url, data))
    puts response
    puts ''
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
end

janusCheck = JanusCheck.new()
#janusCheck.log_in()
#janusCheck.check_log()
#janusCheck.log_out()
#janusCheck.check_admin()
#janusCheck.check_admin_token()
