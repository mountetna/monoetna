# controller.rb

class Controller

  def index(request)

    template = File.read('./server/views/index.html')
    Rack::Response.new(template)
  end

  def upload(request)

    if request.post?

      req = request.POST()
      signature = generate_signature(req)

      if(signature != req['signature'])

        error_message = 'The client is not authorized for this action.'
        # Delete temp file

        #puts request[:params]
        Rack::Response.new({ success: false, error: error_message }.to_json)
      else

        #request[:params][:type]
        # Copy temp file
        Rack::Response.new({ success: true, signature: signature }.to_json)
      end
    else

      Rack::Response.new({ success: false, message: 'Bad request.' }.to_json)
    end
  end

  def persist_solid(request)

  end

  def persist_blob(request)

  end

  # Hash the upload request
  def generate_signature(request)

    ordered_request = Utils.generate_request(request)
    sig = Utils.sign_request(ordered_request, request['algorithm'])
  end

# THIS IS A STUB FOR MAGMA
  def magma_end_point(request)

    if request.post?
      
      params = request.POST()
      time = Time.now.to_i
      req = {
      
        'directory'=> '/ipi/melanoma/',
        'expires'=> 600,
        'algorithm'=> 'MD5',
        'timestamp'=> time,
        'type'=> params['type'],
        'user_email'=> params['user_email'],
        'auth_token'=> params['auth_token'],
        'file_name'=> params['file_name'] + '-ec65f526379d74e8e4dd60a56f37c868',
      }

      ordered_request = Utils.generate_request(req)
      sig = Utils.sign_request(ordered_request, req['algorithm'])

      response = { success: true, request: req, signature: sig }
      Rack::Response.new(response.to_json)
    else

      response = { success: false, message: 'Bad request.' }
      Rack::Response.new(response.to_json)
    end
  end
# THIS IS A STUB FOR MAGMA

end