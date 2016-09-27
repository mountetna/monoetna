# upload_controller.rb
# This controller handles the upload cycle.

class UploadController < Controller

  def initialize_upload()

    Rack::Response.new({ success: false }.to_json())
  end 

  def upload_blob()
    
    Rack::Response.new({ success: false }.to_json())
  end

  def pause_upload()

    Rack::Response.new({ success: false }.to_json())
  end

  def stop_upload()

    Rack::Response.new({ success: false }.to_json())
  end

  def query_upload()

    Rack::Response.new({ success: false }.to_json())
  end
end