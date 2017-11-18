# client_controller.rb
# This controller serves the client pages and code. 

class ClientController

  def initialize(request, action)

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