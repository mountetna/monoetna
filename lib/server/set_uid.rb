class Metis
  class SetUid
    def initialize(app)
      @app = app
    end

    def call(env)
      Metis.instance.tap do |metis|
        @request = Rack::Request.new(env)

        status, headers, body = @app.call(env)

        if (@request.cookies[metis.config(:metis_uid_name)])
          return [ status, headers, body ]
        end

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
  end
end
