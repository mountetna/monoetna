# status_controller.rb
# Provides the status of a file.

class StatusController < Controller

  def status()

    Rack::Response.new({ success: false }.to_json())
  end
end