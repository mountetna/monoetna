# client_controller.rb
# This controller serves the client pages and code. 

class ClientController

  # The redis_service is not required here but it comes along for the ride since
  # we have to implement three args in the other controllers.
  def initialize(redis_service, request, action)

    @redis_service = redis_service
    @request = request
    @params = request.POST()
    @action = action
  end

  def run()  

    return send(@action)
  end

  def index()

    return File.read('./server/views/index.html')
  end

  def user()

    return File.read('./server/views/user.html')
  end
end