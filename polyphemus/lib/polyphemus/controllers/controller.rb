require 'action_controller'

class Polyphemus
  class Controller < Etna::Controller
    def initialize(request, action=nil)
      super
    end
  end

  class StreamingController < Controller

    private

    def can_hijack?
      @request.env['rack.hijack']
    end

    def send_headers(stream)
      headers = [
        "HTTP/1.1 200 OK",
        "Content-Type: text/event-stream"
      ]
      stream.write(headers.map { |header| header + "\r\n" }.join)
      stream.write("\r\n")
      stream.flush
    rescue
      stream.close
      raise
    end
  end
end
