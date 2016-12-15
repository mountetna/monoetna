# admin_controller.rb
# This controller serves as an intermediary to Janus logging services.

class AdminController < Controller

  def initialize(redis_service, request, action)

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

      return send_server_error()
    end 

    if !@user_info['success']

      return send_bad_request()
    end

    # Check that the user has permissions.
    @user_info = @user_info['user_info']
    if !@user_info.key?('permissions')

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

  def add_project()

  end

  def edit_project()

  end

  def delete_project()

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

  def add_permission()

  end

  def edit_projects()

  end

  def delete_permission()

  end

  def make_request(url, data)

    uri = URI.parse(url)
    http = Net::HTTP.new(uri.host, uri.port)
    request = Net::HTTP::Post.new(uri.request_uri)
    request.set_form_data(data)

    begin

      response = http.request(request)
      return Rack::Response.new(response.body)

    rescue Timeout::Error, 
           Errno::EINVAL, 
           Errno::ECONNRESET, 
           EOFError, 
           Net::HTTPBadResponse, 
           Net::HTTPHeaderSyntaxError, 
           Net::ProtocolError => error

      return send_server_error()
    end
  end
end