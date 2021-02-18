require_relative 'user'
require_relative 'hmac'

module Etna
  class Auth
    def initialize(app)
      @app = app
    end

    def call(env)
      request = Rack::Request.new(env)

      # There are three ways to authenticate.
      # Either the route does not require auth,
      # you have an hmac or you have a valid token.
      # Both of these will not validate individual
      # permissions; this is up to the controller
      if [ approve_noauth(request), approve_hmac(request), approve_user(request) ].all?{|approved| !approved}
        return fail_or_redirect(request)
      end

      @app.call(request.env)
    end

    private

    def application
      @application ||= Etna::Application.instance
    end

    def server
      @server ||= application.class.const_get(:Server)
    end

    def params(request)
      request.env['rack.request.params']
    end

    def auth(request, type)
      (etna_param(request, :authorization, nil) || '')[/\A#{type.capitalize} (.*)\z/,1]
    end

    def self.etna_url_param(item)
      :"X-Etna-#{item.to_s.split(/_/).map(&:capitalize).join('-')}"
    end

    def etna_param(request, item, fill='X_ETNA_')
      # This comes either from header variable
      # HTTP_X_SOME_NAME or parameter X-Etna-Some-Name
      #
      # We prefer the param so we can use the header elsewhere

      params(request)[Etna::Auth.etna_url_param(item)] || request.env["HTTP_#{fill}#{item.upcase}"]
    end

    # If the application asks for a redirect for unauthorized users
    def fail_or_redirect(request, msg = 'You are unauthorized')
      return [ 401, { 'Content-Type' => 'text/html' }, [msg] ] unless application.config(:auth_redirect)

      uri = URI(
        application.config(:auth_redirect).chomp('/') + '/login'
      )
      uri.query = URI.encode_www_form(refer: request.url)
      return [ 302, { 'Location' => uri.to_s }, [] ]
    end

    def approve_noauth(request)
      route = server.find_route(request)

      return route && route.noauth?
    end

    def janus_approved?(token)
      return false unless application.config(:janus) && application.config(:janus)[:host]

      janus_client = Etna::Clients::Janus.new(
        token: token,
        ignore_ssl: EtnaApp.instance.config(:ignore_ssl),
        **EtnaApp.instance.config(:janus, environment) || {}
      )
    end

    def approve_user(request)
      token = request.cookies[application.config(:token_name)] || auth(request, :etna)

      return false unless token

      begin
        payload, header = application.sign.jwt_decode(token)

        return false if payload["verify"] && !janus_approved?(token)
        return request.env['etna.user'] = Etna::User.new(payload.map{|k,v| [k.to_sym, v]}.to_h, token)
      rescue
        # bail out if anything goes wrong
        return false
      end
    end

    # Names and order of the fields to be signed.
    def approve_hmac(request)
      hmac_signature = etna_param(request, :signature)

      return false unless hmac_signature

      return false unless headers = etna_param(request, :headers)

      headers = headers.split(/,/).map do |header|
        [ header.to_sym, params(request)[header.to_sym] || etna_param(request, header) ]
      end.to_h

      # Now expect the standard headers
      hmac_params = {
        method: request.request_method,
        host: request.host,
        path: request.path,

        expiration: etna_param(request, :expiration),
        id: etna_param(request, :id),
        nonce: etna_param(request, :nonce),
        headers: headers,
        test_signature: hmac_signature
      }

      begin
        hmac = Etna::Hmac.new(application, hmac_params)
      rescue Exception => e
        return false
      end

      request.env['etna.hmac'] = hmac

      return nil unless hmac.valid?

      # success! set the hmac header params as regular params
      params(request).update(headers)

      return true
    end
  end
end
