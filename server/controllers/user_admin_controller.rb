class UserAdminController < BasicController

  def run()

    m = __method__

    # Check for the correct parameters.
    if !@params.key?('token') then raise_err(:BAD_REQ, 0, m) end

    # Check that the user is an admin (only admins allowed to use this class)
    @data = { :token=> @params['token'], :app_key=> Secrets::APP_KEY }
    uri = '/check-admin-token'
    if !admin_user?(uri, @data) then raise_err(:BAD_REQ, 1, m) end

    # The data being sent back to the client should already be in a JSON format.
    return send(@action)
  end

  def get_users()

    make_request(Conf::JANUS_ADDR+'/get-users', @data)
  end

  def get_projects()

    make_request(Conf::JANUS_ADDR+'/get-projects', @data)
  end

  def get_permissions()

    make_request(Conf::JANUS_ADDR+'/get-permissions', @data)
  end
end