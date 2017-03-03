# This should only server a single user status page.
class ClientController

  def initialize(request, action, logger)

    @request = request
    @params = request.POST()
    @action = action
    @logger = logger
  end

  def run()

    return send(@action)
  end

  def index()

    return File.read('./server/views/index.html')
  end

  def user_admin()

    return File.read('./server/views/user_admin.html')
  end

    def network_utils()

    return File.read('./server/views/network_utils.html')
  end
end