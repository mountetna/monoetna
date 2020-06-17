class UserAdminController < BasicController

  def run()

    # Check for the correct parameters.
    raise_err(:BAD_REQ, 0, __method__) if !@params.key?('token')

    # Check that the user is an admin (only admins allowed to use this class)
    @data = { :token=> @params['token'], :app_key=> Secrets::APP_KEY }
    uri = '/check-admin-token'
    raise_err(:BAD_REQ, 1, __method__) if !admin_user?(uri, @data)

    # The data being sent back to the client should already be in a JSON format.
    return send(@action)
  end

  def get_users()

    make_request(Secrets::JANUS_ADDR+'/get-users', @data)
  end

  def get_projects()

    make_request(Secrets::JANUS_ADDR+'/get-projects', @data)
  end

  def get_permissions()

    make_request(Secrets::JANUS_ADDR+'/get-permissions', @data)
  end

  def upload_permissions()

    raise_err(:BAD_REQ, 0, m) if !@params.key?('permissions')
    @data['permissions'] = @params['permissions']
    make_request(Secrets::JANUS_ADDR+'/upload-permissions', @data)
  end

  def remove_permissions()

    raise_err(:BAD_REQ, 0, m) if !@params.key?('permissions')
    @data['permissions'] = @params['permissions']
    make_request(Secrets::JANUS_ADDR+'/remove-permissions', @data)
  end

  def logout_all()

    make_request(Secrets::JANUS_ADDR+'/logout-all', @data)
  end
end