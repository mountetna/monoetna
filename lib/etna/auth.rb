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
      @params = @request.env['rack.request.params']

      # There are two ways to authenticate.
      # Either you have an hmac or you have
      # a valid token. Both of these will
      # not validate individual permissions;
      # this is up to the controller
      return fail_or_redirect unless create_user || approve_hmac

      @app.call(env)
    end

    private

    def auth(type)
      hmac_param(:authorization,nil)
      (@request.env['HTTP_AUTHORIZATION'] || hmac_param(:authorization) || '')[/\A#{type.capitalize} (.*)\z/,1]
    end

    def create_user
      token = auth(:basic)
      return false unless token

      begin
        payload, header = application.sign.jwt_decode(token)
        return @request.env['etna.user'] = Etna::User.new(payload.map{|k,v| [k.to_sym, v]}.to_h)
      rescue JWT::ExpiredSignature
        return false
      rescue JWT::VerificationError
        return false
      rescue ArgumentError => e
        return false
      end
    end

    def hmac_param(item, fill='X_ETNA_')
      # This comes either from header variable
      # HTTP_X_SOME_NAME or parameter X-Etna-Some-Name

      @request.env["HTTP_#{fill}#{item.upcase}"] || @params[:"X-Etna-#{item.to_s.split(/_/).map(&:capitalize).join('-')}"]
    end

    # Names and order of the fields to be signed.
    def approve_hmac
      hmac_signature = auth(:hmac)

      return false unless hmac_signature

      return false unless headers = hmac_param(:headers)

      headers = headers.split(/,/).map do |header|
        [ header.to_sym, hmac_param(header) ]
      end.to_h

      # Now expect the standard headers
      params = {
        method: @request.request_method,
        host: @request.host,
        path: @request.path,

        timestamp: hmac_param(:timestamp),
        id: hmac_param(:id),
        nonce: hmac_param(:nonce),
        headers: headers
      }

      hmac = Etna::Hmac.new(
        application.sign,
        params
      )

      return false unless hmac && hmac.valid_signature?(hmac_signature)

      # success! set the hmac header params as regular params
      @params.update(headers)

      return @request.env['etna.hmac'] = true
    end

    # If the application asks for a redirect for unauthorized users
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
  end
end
