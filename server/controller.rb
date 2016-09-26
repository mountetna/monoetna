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

        send_bad_request()
      end

      # A hash of the request for HMAC.
      signature = generate_signature(req)

      # Check that the request signature is valid.
      if signature != req['signature']

        send_not_authorized()
      else

        start_upload_sequence(request, signature)
      end
    else

      send_bad_request()
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
    
      send_file_completed(request)
    elsif File.file?(partial_file_name)
    
      resume_file_upload(request, signature)
    else
    
      create_file_status(request, signature)
      create_partial_file(partial_file_name)
      send_active_status(request.POST() , signature, 0)
    end
  end
  
  def resume_file_upload(request, signature)

    # Get the partial file name and bytes.
    full_path = generate_file_path(request)
    partial_file_name = full_path  +'.part'
    byte_count = File.size(partial_file_name)

    # Check that there is a file listing in Redis and get the file status.
    status_key = generate_status_key(request)
    check_status_entry(stat_key)
    file_status = JSON.parse(@redis_service.retrive_file_status(status_key))

    # Check the file for consistancy.
    check_byte_consistancy(file_status, byte_count)
    file_status['current_byte_position'] = byte_count

    send_active_status(file_status, signature, byte_count)
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

  # Check that there is a status listing in Redis for our partial file.
  def check_status_entry(status_key)

    status_present = @redis.check_if_status_present(status_key)
    if status_present.length == 0

      # log A temp file is present, but there was no listing for it in Redis
      send_server_error()
    end
  end

  def create_file_status(request, signature)

    post_params = request.POST() 
    stat_key = generate_status_key(request)
    redis_status = @redis_service.set_file_status(stat_key, post_params.to_json)

    if redis_status != 'ok'

      send_server_error()
    end
  end

  def create_partial_file(partial_file_name)

    partial_file = File.new(partial_file_name, 'w')
    partial_file.close()
  end

#  def upload_blob(request)
#
#    req = request.POST()
#    signature = generate_signature(req)
#
#    if signature != req['signature']
#
#      send_not_authorized()
#    else
#
#      req = regenerate_request(req) 
#      
#      blob_hash = verify_blob_integrity(req, request)
#
#      # Hash the temp file and see if it matches our 'next_blob_hash'
#      Rack::Response.new({ success: true, request: req, hash:blob_hash}.to_json)
#    end
#  end
#
#  def regenerate_request(req)
#
#    request = {
#      
#      'directory'=> req['directory'],
#      'expires'=> req['expires'],
#      'algorithm'=> req['algorithm'],
#      'timestamp'=> req['timestamp'].to_i,
#      'type'=> req['type'],
#      'user_email'=> req['user_email'],
#      'authorization_token'=> req['authorization_token'],
#      'original_name' => req['original_name'],
#      'file_name'=> req['file_name'],
#      'file_size'=> req['file_size'].to_i,
#
#      'signature'=> req['signature'],
#
#      'current_blob_size'=> req['current_blob_size'].to_i,
#      'current_byte_position'=> req['current_byte_position'].to_i,
#      'next_blob_hash'=> req['next_blob_hash'],
#      'next_blob_size'=> req['next_blob_size'].to_i
#    }
#  end
#
#  def verify_blob_integrity(req, request)
#
#    #temp_file_path = request['blob'][:tempfile].path
#
#    # Extract the directory names from the request.
#    target_dir = req['directory']
#    target_file_name = req['file_name']
#    full_path = Conf::ROOT_DIR + target_dir + target_file_name
#
#    # The name of the status file
#    status_file_name = full_path + '/'+ target_file_name + '.json'
#
#    # Open and parse the status file.
#    status_file = File.open(status_file_name, 'r')
#    status_obj = status_file.read()
#    status_hash = JSON.parse(status_obj)
#
#    current_md5_from_file = status_hash['next_blob_hash']
#
#    temp_file_path = request['blob'][:tempfile].path()
#    md5_from_temp_file = Digest::MD5.hexdigest(File.read(temp_file_path))
#    
#    puts current_md5_from_file
#    puts md5_from_temp_file
#
#    return md5_from_temp_file
#  end

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
    file_status = @redis_service.retrive_file_status(status_key)

    byte_count = File.size(full_path)
    response = {
        
        success: true, 
        request: file_status.to_json(),
        byte_count: byte_count,
        status: 'complete'
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