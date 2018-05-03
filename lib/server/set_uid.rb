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
          expires: Time.now+Metis.instance.config(:token_life)
        )

        return response.finish
      end
    end

    def malformed_uid?(uid)
      uid !~ /\A[a-f0-9]+\z/
    end

    def call(env)
      Metis.instance.tap do |metis|
        @request = Rack::Request.new(env)

        existing_uid = @request.cookies[metis.config(:metis_uid_name)]

        if existing_uid && malformed_uid?(existing_uid)
          return cookie_response('Malformed METIS_UID!', 422, {})
        end

        status, headers, body = @app.call(env)

        return [ status, headers, body ] if existing_uid

        return cookie_response(body, status, headers)
      end
    end
  end
end
