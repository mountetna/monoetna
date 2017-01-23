# upload_controller.rb
# This controller handles the upload cycle.

class UploadController < Controller

  def authorize_upload()

    # Check for the POST params.
    if !@request.post?()

      #puts 'POST params are not present.'
      return send_bad_request()
    end

    # Make sure that all the required parameters are present
    params = @request.POST()
    if !authorization_parameters?(params)

      return send_bad_request()
    end

    # Verify that the auth token is valid.
    user_info = validate_token(params['authorization_token'])
    if !user_info.key?('success')

      return send_server_error()
    end 

    if !user_info['success']

      return send_bad_request()
    end

    # Check that the user has the appropriate permission to upload.
    user_info = user_info['user_info']
    if !user_info.key?('permissions')

      return send_server_error()
    end

    project_name = params['project_name']
    project_role = params['project_role']
    user_permissions = user_info['permissions']

    # The first item in the 'project_ids' variable (index 0) is the 'group_id' 
    # that the project belongs to. The second item in the 'project_ids' variable
    # is the project id itself.
    project_ids = fetch_project_id(project_name, project_role, user_permissions)
    if project_ids == nil

      return send_bad_request()
    end

    if project_ids[0] != params['group_id'].to_i()

      return send_bad_request();
    end

    # Check that this file system is in sync with the auth server.
    # If there is a project/permission set in Janus there should be a
    # corresponding directory in Metis. 

    directory = Conf::ROOT_DIR+'/'+project_ids[0].to_s+'/'+project_ids[1].to_s()
    if !File.directory?(directory)

      return send_server_error()
    end

    # Check if the file already exists. There will be another check from the
    # client when the user selects the file or updates the file name. The check
    # here is only for a bit of saftey. We also check on the file status in
    # redis as well to make sure that there isn't a current record.
    full_path = directory + '/' + params['file_name']
    if File.file?(full_path)

      return send_bad_request()
    end

    @status_key = generate_status_key()
    @file_status = @redis_service.retrieve_file_status(@status_key)
    if @file_status != nil

      return send_bad_request()
    end

    return generate_authorization(params, user_info, project_ids[1], directory)
  end

  def generate_authorization(params, user_info, project_id, directory)

    time = Time::now.to_i
    sig_algo = 'MD5'
    
    # The redis index SHOULD be a unique key/index for an entry in redis
    redis_index = @redis_service.get_new_index()
    old_index = params['redis_index']

    params = {
    
      'directory'=> directory,
      'expires'=> Conf::UPLOAD_EXPIRE,
      'signing_algorithm'=> sig_algo,
      'hashing_algorithm'=> sig_algo,
      'start_timestamp'=> time,
      'authorization_token'=> params['authorization_token'],
      'original_name' => params['original_name'],
      'file_name'=> params['file_name'],
      'file_size'=> params['file_size'],
      'user_email'=> user_info['email'],
      'user_id'=> user_info['user_id'],
      'project_id'=> project_id, 
      'old_index'=> old_index,
      'redis_index'=> redis_index,
      'group_id'=> params['group_id']
    }

    ordered_params = SignService::order_params(params)
    sig = SignService::sign_request(ordered_params,sig_algo)
    params['status'] = 'authorized'

    response = { 

      :success=> true,
      :request=> params,
      :signature=> sig,
      :status=> 'authorized' 
    }

    Rack::Response.new(response.to_json)
  end

  def start_upload()

    generate_common_items()

    if !request_valid?()

      return send_bad_request()
    end

    create_file_status()
    create_partial_file()

    if @file_status == nil

      return send_server_error()
    end

    return send_upload_initiated()
  end

  def create_file_status()

    params = @request.POST()
    params['current_blob_size'] = 0
    params['current_byte_position'] = 0
    @redis_service.set_file_status(@status_key, params.to_json)
    @file_status = @redis_service.retrieve_file_status(@status_key)
  end

  def create_partial_file()

    partial_file = File.new(@partial_file_name, 'w')
    partial_file.close()
  end

  def upload_blob()

    generate_common_items()

    if !request_valid?()

      return send_bad_request()
    end

    if @file_status == nil

      return send_server_error()
    end

    if !blob_integrity_ok?()

      return send_bad_request()
    end

    append_blob()
    update_file_status()

    if upload_complete?()

      make_file_permanent()
      hash_and_set_file_status()
      return send_upload_complete()
    else

      return send_upload_active()
    end
  end

  def append_blob()

    temp_file_name = @request['blob'][:tempfile].path()
    partial_file = File.open(@partial_file_name, 'ab')
    temp_file = File.open(temp_file_name, 'rb')
    partial_file.write(temp_file.read())

    partial_file.close()
    temp_file.close()
  end

  def update_file_status()

    params = @request.POST()
    temp_file_path = @request['blob'][:tempfile].path()

    @file_status['current_byte_position'] = File.size(@partial_file_name)
    @file_status['current_blob_size'] = File.size(temp_file_path)
    @file_status['next_blob_hash'] = params['next_blob_hash']
    @file_status['next_blob_size'] = params['next_blob_size']
    @file_status['status'] = 'active'

    @redis_service.set_file_status(@status_key, @file_status.to_json)
  end

    # Generate commonly used variables that we will reuse in many places
  def generate_common_items()

    @signature = generate_signature()
    @status_key = generate_status_key()
    @file_status = @redis_service.retrieve_file_status(@status_key)
    @full_path = generate_file_path()
    @partial_file_name = @full_path  +'.part'
  end

  # Hash the upload request.
  def generate_signature()

    params = @request.POST()
    ordered_params= SignService::order_params(params)
    sig = SignService::sign_request(ordered_params, params['signing_algorithm'])
  end

  # Extract the directory/file names from the request.
  def generate_file_path()

    params = @request.POST()
    return params['directory'] + '/' + params['file_name']
  end

  # Generate the key used to access the file's metadata in Redis.
  def generate_status_key()

    params = @request.POST()
    status_key = params['redis_index'] + '.'
    status_key = status_key + params['file_name'] + '.'
    status_key = status_key + params['group_id'] + '.'
    status_key = status_key + params['project_id']
  end

  # Check the validity of the upload request.
  def request_valid?()

    if upload_errors?()

      return false
    end

    if !@request.post?()

      #puts 'POST params are not present.'
      return false
    end

    if !SignService::request_parameters_valid?(@request.POST())

      #puts 'POST params are not in the correct format.'
      return false
    end

    if generate_signature() != @request.POST()['signature']

      #puts 'The packet doesn\'t have the correct HMAC signature/hash'
      return false
    end

    # We need to make sure that the system clocks of Metis and Magma are in sync
    start_timestamp = @request.POST()['start_timestamp'].to_i
    expiration = @request.POST()['expires'].to_i
    now = Time::now.to_i
    if now >= (start_timestamp + expiration)

      #puts 'The request is past it\'s expiration time.'
      return false
    end
  
    return true
  end

  def authorization_parameters?(params)

    if !params.key?('authorization_token')

      return false
    end

    if !params.key?('project_name')

      return false
    end

    if !params.key?('project_role')

      return false
    end

    if !params.key?('project_id')

      return false
    end

    if !params.key?('file_name')

      return false
    end

    if !params.key?('redis_index')

      return false
    end

    if !params.key?('user_id')

      return false
    end

    if !params.key?('group_id')

      return false
    end

    return true
  end

  def upload_errors?()

    status = @redis_service.retrieve_file_status(@status_key)

    if File.file?(@full_path) && status == nil

      #puts 'FILE_NO_STATUS'
      return true
    end

    if File.file?(@partial_file_name) && status == nil

      #puts 'TEMP_NO_STATUS'
      return true
    end

    if File.file?(@full_path) && File.file?(@partial_file_name)

      #puts 'TEMP_AND_FILE'
      return true
    end

    if !status.nil?

      if !File.file?(@full_path) && !File.file?(@partial_file_name)

        #puts 'STATUS_NO_TEMP_OR_FILE'
        return true
      end
    end

    return false
  end

  def blob_integrity_ok?()

    # Check the blob hash.
    if !blob_hash_ok?()

      return false
    end

    # Check the blob size.
    temp_file_path = @request['blob'][:tempfile].path()
    if File.size(temp_file_path).to_i != @file_status['next_blob_size'].to_i

      return false
    end

    return true
  end

  def blob_hash_ok?()

    md5_from_status = @file_status['next_blob_hash']
    temp_file_path = @request['blob'][:tempfile].path()
    md5_of_temp_file = Digest::MD5.hexdigest(File.read(temp_file_path))

    if md5_from_status != md5_of_temp_file

      return false 
    else

      return true
    end
  end

  def upload_complete?()

    if File.size(@partial_file_name) == @file_status['file_size'].to_i

      return true
    else

      return false
    end
  end

  def make_file_permanent()

    File.rename(@partial_file_name, @full_path)
  end

  def hash_and_set_file_status()

    @file_status.delete('authorization_token')
    @file_status.delete('current_blob_size')
    @file_status.delete('current_byte_position')
    @file_status.delete('expires')
    @file_status.delete('next_blob_size')
    @file_status.delete('next_blob_hash')
    @file_status.delete('signature')
    @file_status.delete('signing_algorithm')
    @file_status.delete('status')

    @file_status['finish_timestamp'] = Time::now.to_i
    @file_status['file_size'] = File.size(@full_path)
    @file_status['hash'] = Digest::MD5.hexdigest(File.read(@full_path))

    @redis_service.set_file_status(@status_key, @file_status.to_json)
  end

  def pause_upload()

    Rack::Response.new({ :success=> false }.to_json())
  end

  def stop_upload()

    Rack::Response.new({ :success=> false }.to_json())
  end
end