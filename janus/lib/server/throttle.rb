require 'rack/throttle'

class Janus
  class Throttle < Rack::Throttle::Daily
    def allowed?(request)
      if request.path =~ %r!/?time-signature!
        super
      else
        true
      end
    end
  end
end
