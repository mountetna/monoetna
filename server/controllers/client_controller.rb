# This should only server a single user status page.
class ClientController

  def initialize(request, action)

    @request = request
    @params = request.POST()
    @action = action
  end

  def run()

    return send(@action)
  end

  def basic_view()

    return File.read('./server/views/basic_view.html')
  end

  def user_admin()

    return File.read('./server/views/user_admin.html')
  end

  def network_utils()

    return File.read('./server/views/network_utils.html')
  end

  def logged_out()

    return File.read('./server/views/logged_out.html')
  end
end