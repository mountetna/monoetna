# admin_controller.rb
# This controller serves as an intermediary to Janus logging services.

class AdminController < Controller

  def initialize(redis_service, request, action, logger)

    @request = request
    @action = action
    @user_info = nil
    @params = nil
  end

  def run()

    # Check for the POST @params.
    if !@request.post?()

      #puts 'POST @params are not present.'
      return send_bad_request()
    end

    @params = @request.POST()
    if !@params.key?('token')

      return send_bad_request()
    end

    # Verify that the auth token is valid.
    @user_info = validate_token(@params['token'])
    if !@user_info.key?('success')

      @logger.warn('Janus returned an invalid response on a token validation.')
      return send_server_error()
    end 

    if !@user_info['success']

      return send_bad_request()
    end

    if !@user_info.key?('user_info')

      return send_bad_request()
    end

    # Check that the user has permissions.
    @user_info = @user_info['user_info']
    if !@user_info.key?('permissions')

      @logger.warn('There is no permission entry for the requested user.')
      return send_server_error()
    end

    # Check for any admin permissions.
    adminPerms = false
    @user_info['permissions'].each do |permission|

      if permission['role'] == 'administrator'

        adminPerms = true
      end
    end

    if !adminPerms

      return send_bad_request()
    end

    send(@action)
  end

  # Check to see if the user is part of the administration project
  def has_master_perms?()

    masterPerms = false
    @user_info['permissions'].each do |permission|

      if permission['project_name'] == 'administration'

        if permission['project_id'] == 1

          if permission['role'] == 'administrator'

            masterPerms = true
          end
        end
      end
    end
    return masterPerms
  end

  def has_permission_items?()

    has_all_items = true
    @params = @request.POST()

    if !@params.key?('id')           then has_all_items = false end
    if !@params.key?('project_id')   then has_all_items = false end
    if !@params.key?('project_name') then has_all_items = false end
    if !@params.key?('react_key')    then has_all_items = false end
    if !@params.key?('role')         then has_all_items = false end
    if !@params.key?('user_email')   then has_all_items = false end
    if !@params.key?('user_id')      then has_all_items = false end

    return has_all_items
  end

  def get_users()

    if !has_master_perms?() 

      return send_bad_request()
    end

    url = Conf::JANUS_ADDR
    url = url + '/get-users'
    data = { 

      :token=> @params['token'], 
      :app_key=> Conf::APP_KEY 
    }
    response = make_request(url, data)
    return response
  end

  def add_users()

    if !has_master_perms?()

      return send_bad_request()
    end
  end

  def edit_user()

    if !has_master_perms?()

      return send_bad_request()
    end
  end

  def delete_user()

    if !has_master_perms?()

      return send_bad_request()
    end
  end

  def get_projects()

    url = Conf::JANUS_ADDR
    url = url + '/get-projects'
    data = { 

      :token=> @params['token'], 
      :app_key=> Conf::APP_KEY 
    }
    response = make_request(url, data)
    return response
  end

  def get_permissions()

    url = Conf::JANUS_ADDR
    url = url + '/get-permissions'
    data = { 

      :token=> @params['token'], 
      :app_key=> Conf::APP_KEY 
    }
    response = make_request(url, data)
    return response
  end

  def save_permission()

    if !has_permission_items?()

      return send_bad_request()
    end

    url = Conf::JANUS_ADDR
    url = url + '/save-permission'
    data = {

      :token=> @params['token'], 
      :app_key=> Conf::APP_KEY,
      :id=> @params['id'],
      :project_id=> @params['project_id'],
      :project_name=> @params['project_name'],
      :react_key=> @params['react_key'],
      :role=> @params['role'],
      :user_email=> @params['user_email'],
      :user_id=> @params['user_id']
    }

    response = make_request(url, data)
    return response
  end

  def upload_permissions()

    if !@params.key?('permissions')

      puts 'sup'
      return send_bad_request()
    end

    url = Conf::JANUS_ADDR
    url = url + '/upload-permissions'
    data = {

      :token=> @params['token'], 
      :app_key=> Conf::APP_KEY,
      :permissions=> @params['permissions']
    }

    response = make_request(url, data)
    return response
  end

  def make_request(url, data)

    uri = URI.parse(url)
    http = Net::HTTP.new(uri.host, uri.port)
    request = Net::HTTP::Post.new(uri.request_uri)
    request.set_form_data(data)

    begin

      response = http.request(request)
      response_code = response.code.to_i()

      if response_code == 200

        return Rack::Response.new(response.body)
      else

        met = __method__.to_s + ', '+ url +', '+ response_code.to_s
        return send_server_error(6, met)
      end
    rescue Timeout::Error, 
           Errno::EINVAL, 
           Errno::ECONNRESET, 
           EOFError, 
           Net::HTTPBadResponse, 
           Net::HTTPHeaderSyntaxError, 
           Net::ProtocolError => error

      met = __method__.to_s + ', '+ url +', '+ response_code.to_s
      return send_server_error(6, met)
    end
  end
end