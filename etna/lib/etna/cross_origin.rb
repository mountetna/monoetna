module Etna
  class CrossOrigin
    def initialize(app)
      @app = app
    end

    def call(env)
      request = Rack::Request.new(env)

      # Don't filter unless this is a CORS request
      return actual_response(request) unless has_header?(request, :origin)

      # OPTIONS requests are a 'preflight' request
      if request.request_method == 'OPTIONS'
        return preflight_response(request)
      end

      # The actual request following a preflight should
      # set some extra headers
      return postflight_response(request)
    end

    private

    def subdomain(host)
      host.split(/\./)[-2..-1]
    end

    def origin_allowed?(request)
      subdomain(URI.parse(header(request, :origin)).host) == subdomain(Etna::Application.instance.host)
    end

    def actual_response(request)
      @app.call(request.env)
    end

    def postflight_response(request)
      status, headers, body = actual_response(request)

      if origin_allowed?(request)
        headers.update(
          'Access-Control-Allow-Origin' => header(request, :origin),
          'Access-Control-Allow-Credentials' => 'true'
        )
      end

      return [ status, headers, body ]
    end

    def preflight_response(request)
      [
        200,
        origin_allowed?(request)  ?  {
          'Access-Control-Allow-Methods' => 'GET, POST, PUT, DELETE, OPTIONS',
          'Access-Control-Allow-Headers' => header(request, :access_control_request_headers),
          'Access-Control-Allow-Origin' => header(request, :origin),
          'Access-Control-Allow-Credentials' => 'true'
        } : {},
        ['']
      ]
    end

    def header(request, name)
      request.get_header("HTTP_#{name.to_s.upcase}")
    end

    def has_header?(request, name)
      request.has_header?("HTTP_#{name.to_s.upcase}")
    end
  end
end
