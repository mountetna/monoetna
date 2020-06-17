class Janus
  class Controller < Etna::Controller
    VIEW_PATH = File.expand_path('../views', __dir__)

    def initialize(request, action)
      super
      @token = nil
    end

    private

    def success_json(hash = {})
      success(hash.to_json, 'application/json')
    end

    def h(s)
      ERB::Util.html_escape(s)
    end

    # Quick check that the email is in a somewhat valid format.
    def email_valid?(eml)
      eml =~ /\A[\w+\-.]+@[a-z\d\-]+(\.[a-z\d\-]+)*\.[a-z]+\z/
    end
  end
end
