# This controller serves as an intermediary to Janus logging services.
class UserLogController < Metis::Controller
  def run()
    m = __method__

    # Depending on whether we get token or email/pass combo we perform different
    # checks.
    unless @action == 'log_in'
      # Check that a token is present.
      raise_err(:BAD_REQ, 0, m) if !@params.key?('token')
    else
      # Check that the email/pass is present.
      raise_err(:BAD_REQ,0,m) if !@params.key?('email') || !@params.key?('pass')
    end

    # The data being sent back to the client should already be in a JSON format.
    return send(@action)
  end

  def log_in()

      m = __method__
      set_login_data()
      return make_request(Secrets::JANUS_ADDR+'/login', @data)
  end

  def check_log()

    m = __method__
    set_log_data()
    return make_request(Secrets::JANUS_ADDR+'/check', @data)
  end

  def log_out()

    set_log_data()
    return make_request(Secrets::JANUS_ADDR+'/logout', @data)
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
