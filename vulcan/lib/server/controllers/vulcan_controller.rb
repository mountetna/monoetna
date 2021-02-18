require 'etna'

class Vulcan
  class Controller < Etna::Controller
    private

    def storage
      @storage ||= Vulcan::Storage.new
    end

    def redirect_to(path)
      @response.redirect(path,302)
      @response.finish
    end

    def token
      @token ||= @request.cookies[Vulcan.instance.config(:token_name)]
    end
  end
end
