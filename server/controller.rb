# controller.rb

class Controller

  def index(request)

    template = File.read('./server/views/index.html')
    Rack::Response.new(template)
  end

  def upload(request)

    signature = generate_signature(request)

    if(signature != request[:params][:signature])

      error_message = 'The client is not authorized for this action.'
      Rack::Response.new({ success: false, error: error_message }.to_json)
    else

      Rack::Response.new({ success: true, signature: signature }.to_json)
    end
  end

  def persist_solid(request)

  end

  def persist_blob(request)

  end

  # Hash the upload request
  def generate_signature(request)

    ordered_request = Utils.generate_request(request[:params])
    signature = Utils.sign_request(ordered_request, request[:params][:algorithm])
  end

# THIS IS A STUB FOR MAGMA
  def magma_end_point(request)

    params = request[:params]
    time = Time.now.to_i

    request = {

      directory: '/ipi/melanoma/',
      expires: 600,
      algorithm: 'MD5',
      timestamp: time,
      type: params['type'],
      user_email: params['user_email'],
      auth_token: params['auth_token'],
      file_name: params['file_name'] + '-ec65f526379d74e8e4dd60a56f37c868',
    }

    ordered_request = Utils.generate_request(request)
    signature = Utils.sign_request(ordered_request, request[:algorithm])

    Rack::Response.new({ success: true, request: request, signature: signature }.to_json)
  end
# THIS IS A STUB FOR MAGMA

end