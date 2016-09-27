# client_controller.rb
# This controller serves the client pages and code. 

class ClientController < Controller

  def run()  

    send(@action)
  end

  def index()

    template = File.read('./server/views/index.html')
    Rack::Response.new(template)
  end
end