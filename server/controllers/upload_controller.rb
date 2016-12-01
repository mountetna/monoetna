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
    project_id = fetch_project_id(project_name, project_role, user_permissions)
    if project_id == nil

      return send_bad_request()
    end

    # Check that this file system is in sync with the auth server.
    # If there is a project/permission set in Janus there should be a
    # corresponding directory in Metis. 
    prj_nm =  params['project_name']
    directory = Conf::ROOT_DIR + '/' + prj_nm
    if !File.directory?(directory)

      return send_server_error()
    end

    # Check if the file already exists. There will be another check from the
    # client when the user selects the file or updates the file name. The check
    # here is only for a bit of saftey. We also check on the file status in
    # redis as well to make sure that there isn't a current record.
    full_path = Conf::ROOT_DIR + '/' + prj_nm + '/' + params['file_name']
    if File.file?(full_path)

      return send_bad_request()
    end

    @status_key = generate_status_key()
    if file_status_ok?

      return send_bad_request()
    end

    return generate_authorization(params, user_info, project_id)
  end

  def start_upload()

    generate_common_items()

    if !request_valid?()

      return send_bad_request()
    end

    create_file_status()
    create_partial_file()

    if !file_status_ok?()

      return send_server_error()
    end

    return send_upload_active()
  end

  def create_file_status()

    params = @request.POST()
    params['current_blob_size'] = 0
    params['current_byte_position'] = 0
    @redis_service.set_file_status(@status_key, params.to_json)
  end

  def create_partial_file()

    create_project_directory()
    partial_file = File.new(@partial_file_name, 'w')
    partial_file.close()
  end

  def create_project_directory()

    req = @request.POST()
    project_id = req['project_id']
  end

  def upload_blob()

    generate_common_items()

    if !request_valid?()

      return send_bad_request()
    end

    if !file_status_ok?()

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

    @redis_service.set_file_status(@status_key, @file_status.to_json)
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