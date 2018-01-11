require_relative 'user'
require_relative 'hmac'

module Etna
  class Auth
    def initialize(app)
      raise "Auth must be the last middleware" unless app.is_a?(Etna::Server)
      @app = app
    end

    def call(env)
      @request = Rack::Request.new(env)

      # There are two ways to authenticate.
      # Either you have an hmac or you have
      # a valid token. Both of these will
      # not validate individual permissions;
      # this is up to the controller
      return fail_or_redirect unless create_user || approve_hmac

      @app.call(env)
    end

    private

    def approve_hmac
      false
    end

    def auth(type)
      (@request.env['HTTP_AUTHORIZATION'] || '')[/\A#{type.capitalize} (.*)\z/,1]
    end

    def create_user
      token = auth(:basic)
      return false unless token

      begin
        payload, header = application.sign.jwt_decode(token)
        @request.env['etna.user'] = Etna::User.new(payload.map{|k,v| [k.to_sym, v]}.to_h)
      rescue JWT::ExpiredSignature
        return false
      rescue JWT::VerificationError
        return false
      rescue ArgumentError => e
        return false
      end
    end

    def fail_or_redirect(msg = 'You are unauthorized')
      return [ 401,{},[msg] ] unless application.config(:auth_redirect)

      uri = URI(
        application.config(:auth_redirect).chomp('/') + '/login'
      )
      uri.query = URI.encode_www_form(refer: @request.url)
      return [ 302, { 'Location' => uri.to_s }, [] ]
    end

    def application
      @application ||= Etna::Application.find(@app.class)
    end

    # Names and order of the fields to be signed.
    def approve_hmac
      hmac_header = auth(:hmac)
      return false unless hmac_header

      hmac, signature = Etna::Hmac.from_request(
        @request,
        application.sign
      )

      return false unless hmac.valid_signature?(signature)
    end
  end
end
