# magma_controller.rb
# This is just a temporary stub for magma and will go away soon

class MagmaController < Controller

  def run()  

    send(@action)
  end

  # THIS IS A STUB FOR MAGMA
  def magma_end_point()

    if @request.post?()
      
      params = @request.POST()
      time = Time::now.to_i
      params = {
      
        'directory'=> '/ipi/bladder/',
        'expires'=> 86400,# set at one day for no good reason outside of testing
        'signing_algorithm'=> 'MD5',
        'hashing_algorithm'=> 'MD5',
        'start_timestamp'=> time,
        'authorization_token'=> params['authorization_token'],
        'original_name' => params['original_name'],
        'file_name'=> 'IPI_ABC_XYZ.fcs',
        'file_size'=> params['file_size'].to_i,
        'user_id'=> 12345,
        'group_id'=> 42
      }

      ordered_params= SignService::order_params(params)
      sig = SignService::sign_request(ordered_params, params['signing_algorithm'])

      response = { success: true, request: params, signature: sig }
      Rack::Response.new(response.to_json)
    else

      response = { success: false, message: 'Bad request.' }
      Rack::Response.new(response.to_json)
    end
  end
# THIS IS A STUB FOR MAGMA
end