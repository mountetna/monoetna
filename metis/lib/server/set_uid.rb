class Metis
  class SetUid
    def initialize(app)
      @app = app
    end

    def cookie_response(body, status, headers)
      Metis.instance.tap do |metis|
        response = Rack::Response.new(body, status, headers)
        response.set_cookie(
          metis.config(:metis_uid_name),
          value: metis.sign.uid,
          path: '/',
          domain: metis.config(:token_domain),
          secure: true,
          same_site: :strict,
          expires: Time.now+Metis.instance.config(:token_life)
        )

        return response.finish
      end
    end

    def valid_uid?(uid)
      uid =~ /\A[a-f0-9]{32}\z/
    end

    def call(env)
      Metis.instance.tap do |metis|
        status, headers, body = @app.call(env)

        request = Rack::Request.new(env)

        existing_uid = request.cookies[metis.config(:metis_uid_name)]

        return [ status, headers, body ] if valid_uid?(existing_uid)
        return [ status, headers, body ] if env['rack.hijack_io']

        return cookie_response(body, status, headers)
      end
    end
  end
end
