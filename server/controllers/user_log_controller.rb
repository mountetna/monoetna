# This controller serves as an intermediary to Janus logging services.
class UserLogController < BasicController

  def run()

    # Depending on whether we get token or email/pass combo we perform different
    # checks.
    unless @action == 'log_in'

      # Check that a token is present.
      if !@params.key?('token') then raise_err(:BAD_REQ, 0, __method__) end
    else

      # Check that the email/pass is present.
      if !@params.key?('email') || !@params.key?('pass')

        raise_err(:BAD_REQ, 0, __method__)
      end
    end

    # The data being sent back to the client should already be in a JSON format.
    return send(@action)
  end

  def log_in()

      m = __method__
      set_login_data()
      return make_request(Conf::JANUS_ADDR+'/login', @data)
  end

  def check_log()

    m = __method__
    set_log_data()
    return make_request(Conf::JANUS_ADDR+'/check', @data)
  end

  def log_out()

    set_log_data()
    return make_request(Conf::JANUS_ADDR+'/logout', @data)
  end

  def set_login_data()

    @data = { 

      :email=> @params['email'], 
      :pass=> @params['pass'], 
      :app_key=> Secrets::APP_KEY 
    }
  end

  def set_log_data()

    @data = {

      :token=> @params['token'],
      :app_key=> Secrets::APP_KEY
    }
  end
end