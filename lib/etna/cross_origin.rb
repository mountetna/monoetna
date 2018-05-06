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
      status, headers, body = @app.call(@request.env)

      Rack::Response.new(body, status, headers)
    end

    def postflight_response
      response = actual_response

      if origin_allowed?
        add_headers(
          response,
          access_control_allow_origin: header(:origin),
          access_control_allow_credentials: 'true'
        )
      end

      return response.finish
    end

    def preflight_response
      response = Rack::Response.new('', 200, {})

      if origin_allowed?
        add_headers(
          response,
          access_control_allow_methods: 'GET, POST, PUT, DELETE, OPTIONS',
          access_control_allow_headers: header(:access_control_request_headers),
          access_control_allow_origin: header(:origin),
          access_control_allow_credentials: 'true'
        )
      end

      return response.finish
    end

    def add_headers(response, headers)
      headers.each do |name, value|
        response.set_header(
          name.to_s.split(/_/).map(&:capitalize).join('-'),
          value
        )
      end
    end

    def header(name)
      @request.get_header("HTTP_#{name.to_s.upcase}")
    end

    def has_header?(name)
      @request.has_header?("HTTP_#{name.to_s.upcase}")
    end
  end
end
