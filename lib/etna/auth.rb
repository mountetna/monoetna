require_relative 'user'

module Etna
  class Auth
    def initialize(app)
      raise "Auth must be the last middleware" unless app.is_a?(Etna::Server)
      @app = app
    end

    def fail_or_redirect(msg)
      return [401,{},msg] unless application.config(:auth_redirect)

      uri = URI(
        application.config(:auth_redirect).chomp('/') + '/login'
      )
      uri.query = URI.encode_www_form(refer: @request.url)
      return [ 302, { 'Location' => uri.to_s }, '' ]
    end

    def application
      @application ||= Etna::Application.find(@app.class)
    end

    def call(env)
      @request = Rack::Request.new(env)

      auth = env['HTTP_AUTHORIZATION']

      return fail_or_redirect('Authorization header missing') unless auth

      token = auth.match(/\ABasic (?<token>.*)\z/)[:token]

      return fail_or_redirect('Authorization header malformed') unless token

      application = Etna::Application.find(@app.class)

      begin
        payload, header = application.sign.jwt_decode(token)
        env['etna.user'] = Etna::User.new(payload.map{|k,v| [k.to_sym, v]}.to_h)
      rescue JWT::ExpiredSignature
        return fail_or_redirect('Token signature is expired')
      rescue JWT::VerificationError
        return fail_or_redirect('Token verification failed.')
      rescue ArgumentError
        return fail_or_redirect('No email id set in token.')
      end

      @app.call(env)
    end
  end
end
