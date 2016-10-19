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

    if !request_valid?()

      return send_bad_request()
    end
  
    generate_common_items()

    if upload_errors?()

      return send_server_error()
    end

    return send(@action)
  end

  # Check the validity of the request.
  def request_valid?()

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

  def upload_errors?()

    status = @redis_service.retrive_file_status(@status_key)

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

  def file_status_ok?()

    @file_status = @redis_service.retrive_file_status(@status_key)
    if @file_status == nil

      return false
    else

      @file_status = JSON.parse(@file_status)
      return true
    end
  end

  # Generate commonly used variables that we will reuse in many places
  def generate_common_items()

    @signature = generate_signature()
    @status_key = generate_status_key()
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

    req = @request.POST()
    full_path = Conf::ROOT_DIR + req['directory'] + req['file_name']
  end

  # Generate the key used to access the file's metadata in Redis.
  def generate_status_key()

    req = @request.POST()
    status_key = req['redis_index'] +'.'
    status_key = status_key + req['file_name'] +'.'
    status_key = status_key + req['group_id'] +'.'
    status_key = status_key + req['user_id']
  end

  def send_bad_request()

    Rack::Response.new({ success: false, error: 'Bad request.' }.to_json())
  end

  def send_server_error()

    error_message = 'There was a server error.'
    Rack::Response.new({ success: false, error: error_message }.to_json())
  end

  def send_upload_active()

    response = {

      success: true,
      request: @file_status,
      signature: @signature,
      byte_count: File.size(@partial_file_name),
      status: 'active'
    }
    return Rack::Response.new(response.to_json)
  end

  def send_upload_complete()

    response = {

      success: true,
      result: @file_status,
      status: 'complete'
    }
    return Rack::Response.new(response.to_json())
  end
end