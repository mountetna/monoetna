# controller.rb

class Controller

  def index(request)

    template = File.read('./server/views/index.html')
    Rack::Response.new(template)
  end

  def upload(request)

    if request.post?()

      req = request.POST()
      signature = generate_signature(req)

      if signature != req['signature']

        remove_temp_file(req)
        error_message = 'The client is not authorized for this action.'
        Rack::Response.new({ success: false, error: error_message }.to_json)
      else

        persist_solid(req)
        Rack::Response.new({ success: true, signature: signature }.to_json)
      end
    else

      Rack::Response.new({ success: false, error: 'Bad request.' }.to_json)
    end
  end

  def remove_temp_file(request)

    temp_file_path = request['upload_file'][:tempfile].path
    if File::exists?(temp_file_path)

      File::delete(temp_file_path)
    end
  end

  def persist_solid(request)

    # Extract the directory names from the request.
    temp_file_path = request['upload_file'][:tempfile].path
    target_dir = request['directory']
    target_file_name = request['file_name']

    # Create the target directory if it does not exist and move the file.
    dir = Conf::ROOT_DIR + target_dir
    create_directory(dir, target_dir)
    FileUtils::mv(temp_file_path, Conf::ROOT_DIR+target_dir+target_file_name)
  end

  def persist_blob(request)

  end

  def create_directory(dir, target_dir)

    if !File::directory?(dir)

      FileUtils::mkdir_p(Conf::ROOT_DIR + directory)
    end
  end

  # Hash the upload request
  def generate_signature(request)

    ordered_request = Utils::generate_request(request)
    sig = Utils::sign_request(ordered_request, request['algorithm'])
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
        'type'=> params['type'],
        'user_email'=> params['user_email'],
        'auth_token'=> params['auth_token'],
        'file_name'=> params['file_name'],
        'file_size'=> params['file_size'].to_i
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