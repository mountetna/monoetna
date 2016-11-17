# magma_controller.rb
# This is just a temporary stub for magma and will go away soon

class MagmaController < Controller

  def run()  

    send(@action)
  end

  def check_auth(auth_token)

    uri = URI('http://janus-stage.ucsf.edu/check')
    request = Net::HTTP::Post.new(uri)
    request.set_form_data({ :token=> auth_token, :app_key=> Conf::APP_KEY })

    response = Net::HTTP.start(uri.hostname, uri.port) do |http|

      http.request(req)
    end

    case response
    when Net::HTTPSuccess, Net::HTTPRedirection

      puts response.value
    else

      puts response.value
    end
  end

  # THIS IS A STUB FOR MAGMA
  def magma_end_point()

    if @request.post?()
      
      params = @request.POST()

      auth_token = params['authorization_token']
      # Do a auth cycle here with the user's auth token.

      check_auth(auth_token)

      # The redis index SHOULD be a unique key/index for an entry in redis
      redis_index = @redis_service.get_new_index()
      old_index = params['redis_index']

      time = Time::now.to_i
      params = {
      
        'directory'=> '/ipi/bladder/',
        'expires'=> 86400,# set at one day for no good reason outside of testing
        'signing_algorithm'=> 'MD5',
        'hashing_algorithm'=> 'MD5',
        'start_timestamp'=> time,
        'authorization_token'=> params['authorization_token'],
        'original_name' => params['original_name'],
        'file_name'=> params['original_name'],
        'file_size'=> params['file_size'].to_i,
        'user_email'=> params['user_email'],
        'user_id'=> 12345,
        'group_id'=> 42, 
        'old_index'=> old_index,
        'redis_index'=> redis_index
      }

      ordered_params= SignService::order_params(params)
      sig_algo = params['signing_algorithm']
      sig = SignService::sign_request(ordered_params,sig_algo)

      response = { success: true, request: params, signature: sig }
      Rack::Response.new(response.to_json)
    else

      response = { success: false, message: 'Bad request.' }
      Rack::Response.new(response.to_json)
    end
  end
# THIS IS A STUB FOR MAGMA
end