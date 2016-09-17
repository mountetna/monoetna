# controller.rb

class Controller

  def index(request)

    template = File.read('./server/views/index.html')
    Rack::Response.new(template)
  end

#  def upload(request)
#
#    if request.post?()
#
#      req = request.POST()
#      signature = generate_signature(req)
#
#      if signature != req['signature']
#
#        remove_temp_file(req)
#      else
#
#        persist_solid(req)
#        Rack::Response.new({ success: true, signature: signature }.to_json)
#      end
#    else
#
#      send_bad_request()
#    end
#  end

  def initialize_upload(request)

    #1 check the sig/hash
    #2 check for file or directory
      #2.1 make file directory and save status card
    #3 send status back to the client

    if request.post?()

      req = request.POST()
      signature = generate_signature(req)

      if signature != req['signature']

        send_not_authorized()
      else

        # Extract the directory names from the request.
        target_dir = request['directory']
        target_file_name = request['file_name']
        complete_path = Conf::ROOT_DIR + target_dir + target_file_name

        if File.file?(complete_path)

          send_file_complete(req, signature, complete_path)
        else

          # The name of the status file
          status_file_name = complete_path + '/'+ target_file_name + '.json'

          if File.directory?(complete_path)

            status_file = File.open(status_file_name, 'r')
            status_obj = status_file.read()
            status_hash = JSON.parse(status_obj)

            # Check for partial file
            # If no partial file exsists then set byte count to 0
            # Otherwise byte_count equals the partial file size
            partial_file_name = complete_path + '/'+ target_file_name + '.part'
            byte_count = 0
            if File.file?(partial_file_name)

              byte_count = File.size(partial_file_name)
            end

            response = {

              success: true,
              request: req.to_json(),
              signature: signature,
              byte_count: byte_count,
              status: 'active'
            }

            Rack::Response.new(response.to_json)
          else

            # Create the blobs' parent directory if it does not exists.
            create_directory(complete_path)

            #create the first json file
            status_file = File.open(status_file_name, 'w')
            status_file.puts(req.to_json)
            status_file.close()

            response = {

              success: true,
              request: req.to_json(),
              signature: signature,
              byte_count: 0,
              status: 'initialized'
            }

            Rack::Response.new(response.to_json)
          end

          #Rack::Response.new({ exsits: false, dir: complete_path }.to_json)
        end

        #dir = Conf::ROOT_DIR + target_dir + target_file_name
        #Rack::Response.new({ success: true, dir: dir }.to_json)

        # Check to see if the file already exsists (is complete)
#       if File.exist?(Conf::ROOT_DIR + target_dir + target_file_name)

#         #send_file_complete(request, signature)
#       else

#         Rack::Response.new({ success: true, signature: signature }.to_json)
#       end

        # Create the target directory if it does not exist.
        #dir = Conf::ROOT_DIR + target_dir
        #create_directory(dir, target_dir)

        # Check if the file OR file directory exists.

        # If the file directory exists, get the status card and return.
        #if !File::directory?(dir)

        #end

        # If the file exsists, return the status that the file is already
        # present.

        # If the file directory does not exist, then create the temp file dir
        # and create the first status card.
      end
    else

      send_bad_request()
    end
  end

  def upload_blob(request)

    req = request.POST()
    signature = generate_signature(req)

    if signature != req['signature']

      send_not_authorized()
    else

      Rack::Response.new({ success: true, request: req }.to_json)
    end
  end

#  def status(request)
#
#    #check for file or directory
#    #send status back to the client
#  end

  def remove_temp_file(request)

    temp_file_path = request['upload_file'][:tempfile].path()
    if File::exists?(temp_file_path)

      File::delete(temp_file_path)
    end
  end

# def persist_solid(request)

#   # Extract the directory names from the request.
#   temp_file_path = request['upload_file'][:tempfile].path
#   target_dir = request['directory']
#   target_file_name = request['file_name']

#   # Create the target directory if it does not exist and move the file.
#   dir = Conf::ROOT_DIR + target_dir
#   create_directory(dir, target_dir)
#   FileUtils::mv(temp_file_path, Conf::ROOT_DIR+target_dir+target_file_name)
# end

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

  def send_file_complete(request, signature, complete_path)

    byte_count = File.size(complete_path)
    response = {
        
        success: true, 
        request: request.to_json(),
        signature: signature,
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
        'authorization_token'=> params['authorization_token'],
        'original_name' => params['original_name'],
        'file_name'=> params['original_name'] +'-f8cfc63531b6ed753bd536f6e12d578c',
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