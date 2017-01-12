# controller.rb
# The generic controller that handles validations and common processing tasks.

class Controller

  def initialize(redis_service, request, action)

    @redis_service = redis_service
    @request = request
    @action = action

    # see generate_common_items()
    @signature = nil
    @status_key = nil
    @full_path = nil
    @partial_file_name = nil
    @file_status = nil
  end

  def run()

    return send(@action)
  end

  def validate_token(token)

    url = Conf::JANUS_ADDR
    url = url + '/check'
    data = { :token=> token, :app_key=> Conf::APP_KEY }
    #response = make_request(url, data)

    uri = URI.parse(url)
    http = Net::HTTP.new(uri.host, uri.port)
    request = Net::HTTP::Post.new(uri.request_uri)
    request.set_form_data(data)

    begin

      response = http.request(request)
      response_code = response.code.to_i()

      if response_code == 200

        return JSON.parse(response.body)
      else

        return {}
      end
    rescue Timeout::Error, 
           Errno::EINVAL, 
           Errno::ECONNRESET, 
           EOFError, 
           Net::HTTPBadResponse, 
           Net::HTTPHeaderSyntaxError, 
           Net::ProtocolError => error

      return { 'success'=> false }
    end
  end

  def fetch_project_id(project_name, project_role, user_permissions)

    project_ids = nil
    for permission in user_permissions

      if project_name == permission['project_name']

        if project_role == permission['role']

          if project_role == 'administration' || project_role == 'editor'

            project_ids = [permission['group_id'], permission['project_id']]
            break
          end
        end
      end
    end
    return project_ids
  end

  def file_status_ok?()

    @file_status = @redis_service.retrieve_file_status(@status_key)
    if @file_status == nil

      return false
    else

      return true
    end
  end

  def retrieve_files()

    # Check for the POST params.
    if !@request.post?()

      #puts 'POST params are not present.'
      return send_bad_request()
    end

    params = @request.POST()
    if !params.key?('authorization_token')

      return false
    end

    # Verify that the auth token is valid.
    user_info = validate_token(params['authorization_token'])
    if !user_info.key?('success')

      return send_server_error()
    end 

    if !user_info['success']

      return send_bad_request()
    end

    # Extract the projects that the user has permissions on.
    user_info = user_info['user_info']
    if !user_info.key?('permissions')

      return send_server_error()
    end

    project_ids = []
    user_info['permissions'].each do |permission|

      project_id = permission['project_id'].to_s()
      group_id = permission['group_id'].to_s()
      project_ids.push(group_id+'.'+project_id)
    end

    file_list = pull_file_metadata(project_ids)
    Rack::Response.new({ :success=> true, :file_list=> file_list }.to_json())
  end

  def pull_file_metadata(project_ids)

    file_data = []
    for project_id in project_ids

      keys = @redis_service.retrieve_file_key('*'+project_id)
      if keys != nil

        for key in keys

          file_metadata = @redis_service.retrieve_file_status(key)
          if file_metadata != nil
            
            file_data.push(JSON.parse(file_metadata))
          end
        end
      end
    end
    return file_data
  end

  def send_bad_request()

    Rack::Response.new({ :success=> false, :error=> 'Bad request.' }.to_json())
  end

  def send_server_error()

    error_message = 'There was a server error.'
    Rack::Response.new({ :success=> false, :error=> error_message }.to_json())
  end

  def send_upload_active()

    response = {

      :success=> true,
      :request=> @file_status,
      :signature=> @signature,
      :byte_count=> File.size(@partial_file_name),
      :status=> 'active'
    }
    return Rack::Response.new(response.to_json)
  end

  def send_upload_complete()

    response = {

      :success=> true,
      :result=> @file_status,
      :status=> 'complete'
    }
    return Rack::Response.new(response.to_json())
  end
end