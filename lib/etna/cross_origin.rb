module Etna
  class CrossOrigin
    def initialize(app)
      @app = app
    end

    def call(env)
      @request = Rack::Request.new(env)

      # Don't filter unless this is a CORS request
      return actual_response unless has_header?(:origin)

      # OPTIONS requests are a 'preflight' request
      if @request.request_method == 'OPTIONS'
        return preflight_response
      end

      # The actual request following a preflight should
      # set some extra headers
      return postflight_response
    end

    private

    def subdomain(host)
      host.split(/\./)[-2..-1]
    end

    def origin_allowed?
      subdomain(URI.parse(header(:origin)).host) == subdomain(@request.host)
    end

    def actual_response
      @app.call(@request.env)
    end

    def postflight_response
      status, headers, body = actual_response

      if origin_allowed?
        headers.update(
          'Access-Control-Allow-Origin' => header(:origin),
          'Access-Control-Allow-Credentials' => 'true'
        )
      end

      return [ status, headers, body ]
    end

    def preflight_response
      [
        200,
        origin_allowed?  ?  {
          'Access-Control-Allow-Methods' => 'GET, POST, PUT, DELETE, OPTIONS',
          'Access-Control-Allow-Headers' => header(:access_control_request_headers),
          'Access-Control-Allow-Origin' => header(:origin),
          'Access-Control-Allow-Credentials' => 'true'
        } : {},
        ''
      ]
    end

    def header(name)
      @request.get_header("HTTP_#{name.to_s.upcase}")
    end

    def has_header?(name)
      @request.has_header?("HTTP_#{name.to_s.upcase}")
    end
  end
end
