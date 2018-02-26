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
      return fail_or_redirect unless approve_hmac || create_user

      @app.call(env)
    end

    private

    def auth(type)
      (etna_param(:authorization, nil) || '')[/\A#{type.capitalize} (.*)\z/,1]
    end

    def create_user
      token = @request.cookies[application.config(:token_name)] || auth(:etna)

      return false unless token

      begin
        payload, header = application.sign.jwt_decode(token)
        return @request.env['etna.user'] = Etna::User.new(payload.map{|k,v| [k.to_sym, v]}.to_h)
      rescue
        # bail out if anything goes wrong
        return false
      end
    end

    def etna_param(item, fill='X_ETNA_')
      # This comes either from header variable
      # HTTP_X_SOME_NAME or parameter X-Etna-Some-Name
      #
      # We prefer the param so we can use the header elsewhere

      @params[:"X-Etna-#{item.to_s.split(/_/).map(&:capitalize).join('-')}"] || @request.env["HTTP_#{fill}#{item.upcase}"] 
    end

    # Names and order of the fields to be signed.
    def approve_hmac
      hmac_signature = auth(:hmac)

      return false unless hmac_signature

      return false unless headers = etna_param(:headers)

      headers = headers.split(/,/).map do |header|
        [ header.to_sym, etna_param(header) ]
      end.to_h

      # Now expect the standard headers
      params = {
        method: @request.request_method,
        host: @request.host,
        path: @request.path,

        expiration: etna_param(:expiration),
        id: etna_param(:id),
        nonce: etna_param(:nonce),
        headers: headers
      }

      begin
        hmac = Etna::Hmac.new(
          application,
          params
        )
      rescue Exception => e
        return false
      end

      return false unless hmac.valid_signature?(hmac_signature)

      # success! set the hmac header params as regular params
      @params.update(headers)

      return @request.env['etna.hmac'] = true
    end

    # If the application asks for a redirect for unauthorized users
    def fail_or_redirect(msg = 'You are unauthorized')
      return [ 401, { 'Content-Type' => 'text/html' }, [msg] ] unless application.config(:auth_redirect)

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
