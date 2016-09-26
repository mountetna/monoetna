# controller.rb

class Controller

  def initialize(redis_service)

    @redis_service = redis_service
  end

  def index(request)

    template = File.read('./server/views/index.html')
    Rack::Response.new(template)
  end

  # Make sure the client request looks good.
  def initialize_upload(request)

    # Check to see if the POST values are present.
    if request.post?()

      # The POST data from the client.
      req = request.POST() 
      if !Utils::verify_request_parameters(req)

        return send_bad_request()
      end

      # A hash of the request for HMAC.
      signature = generate_signature(req)

      # Check that the request signature is valid.
      if signature != req['signature']

        return send_not_authorized()
      else

        return start_upload_sequence(request, signature)
      end
    else

      return send_bad_request()
    end
  end

  # Check the status of the file that the client is trying to upload.
  def start_upload_sequence(request, signature)

    # Complete path of the file.
    full_path = generate_file_path(request) 

    # Redis key to access the file metadata.
    status_key = generate_status_key(request)

    # Path of a possible partial file.
    partial_file_name = full_path  +'.part'

    # Check to see: if the file to upload already exists; if partial file has 
    # already been started but previously interrupted or stopped; if we need 
    # to start a fresh file status.
    if File.file?(full_path)
    
      return send_file_completed(request)
    elsif File.file?(partial_file_name)
    
      return resume_file_upload(request, signature)
    else
    
      if !create_file_status(request, signature)

        return send_server_error()
      else

        create_partial_file(partial_file_name)
        return send_active_status(request.POST() , signature, 0)
      end
    end
  end
  
  def resume_file_upload(request, signature)

    # Get the partial file name and bytes.
    full_path = generate_file_path(request)
    partial_file_name = full_path  +'.part'
    byte_count = File.size(partial_file_name)

    # Check that there is a file listing in Redis and get the file status.
    status_key = generate_status_key(request)
    
    if !@redis_service.status_present?(status_key)

      # log 'A temp file is present, but there was no listing for it in Redis'
      return send_server_error()
    end

    file_status = JSON.parse(@redis_service.retrive_file_status(status_key))

    if file_status == nil

      # the status return from Redis is nil. That should not happen here.
      return send_server_error()
    end

    # Check the file for consistancy.
    check_byte_consistancy(file_status, byte_count)
    file_status['current_byte_position'] = byte_count

    return send_active_status(file_status, signature, byte_count)
  end

  # Check the diff between file status byte count and the actual file.
  def check_byte_consistancy(file_status, byte_count)

    if !file_status.key?('current_byte_position')

      # log partial_file_name did not have a current bytes listed in it's status
    end

    if byte_count != file_status['current_byte_position']

      # log the partial file size did not match it's listing in it's status
    end
  end

  def create_file_status(request, signature)

    post_params = request.POST() 
    stat_key = generate_status_key(request)
    redis_status = @redis_service.set_file_status(stat_key, post_params.to_json)

    if redis_status != 'ok'

      return true
    else

      return false
    end
  end

  def create_partial_file(partial_file_name)

    partial_file = File.new(partial_file_name, 'w')
    partial_file.close()
  end

  def upload_blob(request)

    req = request.POST()
    signature = generate_signature(req)

    if signature != req['signature']

      send_not_authorized()
    else
      
      # Hash the temp file and see if it matches our 'next_blob_hash'
      if !blob_integrity_ok?(request)

        Rack::Response.new({ success: false }.to_json)
      else

        append_blob(request)

        full_path = generate_file_path(request)
        partial_file_name = full_path  +'.part'

        byte_count = File.size(partial_file_name)

        # Check that there is a file listing in Redis and get the file status.
        status_key = generate_status_key(request)

        file_status = JSON.parse(@redis_service.retrive_file_status(status_key))
        file_status['current_byte_position'] = byte_count
        file_status['current_blob_size'] = req['current_blob_size']
        file_status['next_blob_hash'] = req['next_blob_hash']
        file_status['next_blob_size'] = req['next_blob_size']

        @redis_service.set_file_status(status_key, file_status.to_json)

        if byte_count == file_status['file_size'].to_i

          File.rename(partial_file_name, full_path)
          return send_file_completed(request)
        else
          
          return send_active_status(file_status, signature, byte_count)
        end
      end
    end
  end

  def append_blob(request)

    full_path = generate_file_path(request)
    partial_file_name = full_path  +'.part'
    temp_file_name = request['blob'][:tempfile].path()

    partial_file = File.open(partial_file_name, 'ab')
    temp_file = File.open(temp_file_name, 'rb')

    partial_file.write(temp_file.read())
  end

  def pause_upload(request)

  end

  def stop_upload(request)

  end

  def query_upload(request)

  end

  def blob_integrity_ok?(request)

    status_key = generate_status_key(request)
    
    if !@redis_service.status_present?(status_key)

      # log A temp file is present, but there was no listing for it in Redis
      return send_server_error()
    end

    file_status = JSON.parse(@redis_service.retrive_file_status(status_key))

    if file_status == nil

      # the status return from Redis is nil. That should not happen here.
      return send_server_error()
    end

    md5_from_status = file_status['next_blob_hash']

    temp_file_path = request['blob'][:tempfile].path()
    md5_of_temp_file = Digest::MD5.hexdigest(File.read(temp_file_path))

    if md5_from_status == md5_of_temp_file

      return true
    else

      return false
    end
  end

  def remove_temp_file(request)

    temp_file_path = request['upload_file'][:tempfile].path()
    if File::exists?(temp_file_path)

      File::delete(temp_file_path)
    end
  end

  def create_directory(directory)

    if !File::directory?(directory)

      FileUtils::mkdir_p(directory)
    end
  end

  # Hash the upload request
  def generate_signature(request)

    ordered_request = Utils::generate_request(request)
    sig = Utils::sign_request(ordered_request, request['algorithm'])
  end

  # Extract the directory/file names from the request.
  def generate_file_path(req)

    full_path = Conf::ROOT_DIR + req['directory'] + req['file_name']
  end

  # Generate the key used to access the file's metadata in Redis.
  def generate_status_key(req)

    status_key = req['file_name'] +'.'+ req['user_id'] +'.'+ req['group_id']
  end

  def send_active_status(file_status, signature, byte_count)

    response = {

      success: true,
      request: file_status,
      signature: signature,
      byte_count: byte_count,
      status: 'active'
    }

    Rack::Response.new(response.to_json)
  end

  def send_file_completed(request)

    # Complete path of the file.
    full_path = generate_file_path(request) 

    # Retrieve the file status metadata from Redis
    status_key = generate_status_key(request)
    
    # Check that there is a file listing in Redis and get the file status.
    if !@redis_service.status_present?(status_key)

      # log A temp file is present, but there was no listing for it in Redis
      return send_server_error()
    end

    file_status = @redis_service.retrive_file_status(status_key)

    if file_status == nil

      # the status return from Redis is nil. That should not happen here.
      return send_server_error()
    end

    # Check that the file exists
    if !File::exists?(full_path)

      # Log, we are sending a "completed" response but the file doesn't exsist
      # Maybe it wasn't renamed or was accidentally delete before we got here.
      return send_server_error()
    end

    byte_count = File.size(full_path)
    response = {
        
        success: true, 
        request: file_status,
        byte_count: byte_count,
        status: 'complete',
        signature: file_status['signature']
    }

    Rack::Response.new(response.to_json())
  end

  def send_bad_request()

    Rack::Response.new({ success: false, error: 'Bad request.' }.to_json())
  end

  def send_not_authorized()

    error_message = 'The client is not authorized for this action.'
    Rack::Response.new({ success: false, error: error_message }.to_json())
  end

  def send_server_error()

    error_message = 'There was a server error.'
    Rack::Response.new({ success: false, error: error_message }.to_json())
  end

# THIS IS A STUB FOR MAGMA
  def magma_end_point(request)

    if request.post?
      
      params = request.POST
      time = Time::now.to_i
      req = {
      
        'directory'=> '/ipi/bladder/',
        'expires'=> 600,
        'algorithm'=> 'MD5',
        'timestamp'=> time,
        'user_email'=> params['user_email'],
        'authorization_token'=> params['authorization_token'],
        'original_name' => params['original_name'],
        'file_name'=> 'IPI_ABC_XYZ.fcs',
        'file_size'=> params['file_size'].to_i,
        'user_id'=> 12345,
        'group_id'=> 42
      }

      ordered_request = Utils::generate_request(req)
      sig = Utils::sign_request(ordered_request, req['algorithm'])

      response = { success: true, request: req, signature: sig }
      Rack::Response.new(response.to_json)
    else

      response = { success: false, message: 'Bad request.' }
      Rack::Response.new(response.to_json)
    end
  end
# THIS IS A STUB FOR MAGMA
end