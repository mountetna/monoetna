# client_controller.rb
# This controller serves the client pages and code. 

class ClientController

  # The redis_service is not required here but it comes along for the ride since
  # we have to implement three args in the other controllers.
  def initialize(redis_service, request, action)

    @request = request
    @action = action
  end

  def run()  

    send(@action)
  end

  def index()

    template = File.read('./server/views/index.html')
    Rack::Response.new(template)
  end

  def admin()

    template = File.read('./server/views/admin.html')
    Rack::Response.new(template)
  end

  def user()

    template = File.read('./server/views/user.html')
    Rack::Response.new(template)
  end
end