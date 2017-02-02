# controller.rb
# The generic controller that handles validations and common processing tasks.

class Controller

  def initialize(redis_service, request, action, logger)

    @redis_service = redis_service
    @request = request
    @action = action
    @logger = logger

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

    uri = URI.parse(url)
    http = Net::HTTP.new(uri.host, uri.port)
    request = Net::HTTP::Post.new(uri.request_uri)
    request.set_form_data(data)

    begin

      response = http.request(request)
      response_code = response.code.to_i()

      if response_code == 200

        return user_valid?(JSON.parse(response.body))
      else

        return false
      end
    rescue Timeout::Error, 
           Errno::EINVAL, 
           Errno::ECONNRESET, 
           EOFError, 
           Net::HTTPBadResponse, 
           Net::HTTPHeaderSyntaxError, 
           Net::ProtocolError => error

      return false
    end
  end

  def user_valid?(user_info)

    if !user_info.key?('success')

      return false
    end 

    if !user_info['success']

      return false
    end

    if !user_info.key?('logged')

      return false
    end

    if !user_info['logged']

      return false
    end

    user_info = user_info['user_info']
    if !user_info.key?('permissions')

      return false
    end

    return user_info
  end

  def fetch_project_id(project_name, role, user_permissions)

    project_ids = nil
    for permission in user_permissions

      if project_name == permission['project_name']

        if role == permission['role']

          if role == 'administration' || role == 'editor'

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

      return send_bad_request()
    end

    # Validate the user token.
    user_info = validate_token(params['authorization_token'])
    if !user_info

      return send_bad_request()
    end

    project_ids = []
    user_info['permissions'].each do |permission|

      project_id = permission['project_id'].to_s()
      group_id = permission['group_id'].to_s()
      project_ids.push(group_id+'.'+project_id)
    end

    file_list = pull_file_metadata(project_ids, user_info['user_id'])
    Rack::Response.new({ :success=> true, :file_list=> file_list }.to_json())
  end

  def pull_file_metadata(project_ids, user_id)

    file_data = []
    for project_id in project_ids

      # Pull all the metadata keys that the user has permissions on.
      keys = @redis_service.retrieve_file_key('*'+project_id)
      if keys != nil

        for key in keys

          # Pull the indivitual file status.
          file_metadata = @redis_service.retrieve_file_status(key)

          if file_metadata != nil

            # If the file is complete it will have a 'finished_timestamp'.
            if file_metadata.key?('finish_timestamp')

              file_data.push(file_metadata)
            else

              # If the file is not complete and the 'user_id' matches the user
              # requesting this action, then we know that this is an incomplete
              # upload. We can possibly resume it so we send that back.
              if file_metadata['user_id'].to_s == user_id.to_s

                file_data.push(file_metadata)
              end
            end
          end
        end
      end
    end

    return file_data
  end

  def send_bad_request(id, method)

    ref_id = SecureRandom.hex(4)
    code = Conf::WARNS[id].to_s
    @logger.warn(ref_id.to_s+' - '+code+', '+method.to_s)
    response = {

      :success=> false,
      :error=> 'Bad request.',
      :reference_id=> ref_id,
      :error_code=> id
    }
    Rack::Response.new(response.to_json())
  end

  def send_server_error(id, method)

    ref_id = SecureRandom.hex(4)
    code = Conf::ERRORS[id].to_s
    @logger.error(ref_id.to_s+' - '+code+', '+method.to_s)
    response = { 

      :success=> false,
      :error=> 'There was a server error.',
      :reference_id=> ref_id,
      :error_code=> id
    }
    Rack::Response.new(response.to_json())
  end

  def send_upload_initiated()

    @file_status['status'] = 'initialized'
    response = {

      :success=> true,
      :request=> @file_status,
      :signature=> @signature,
      :byte_count=> File.size(@partial_file_name),
      :status=> 'initialized'
    }

    return Rack::Response.new(response.to_json)
  end

  def send_upload_active()

    @file_status['status'] = 'active'
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

    result = @redis_service.retrieve_file_status(@status_key)
    response = {

      :success=> true,
      :result=> result,
      :status=> 'complete'
    }
    return Rack::Response.new(response.to_json())
  end
end