require_relative 'user'

# This is an authentication layer you can use in testing. It will make an
# Etna::User as usual that your controller can respond to; you can pass
# permissions for this user directly in the Authorization header
module Etna
  class TestAuth
    def initialize(app)
      @app = app
    end

    def application
      @application ||= Etna::Application.find(@app.class)
    end

    def failure(status, msg)
      [ status, {}, [ msg] ]
    end

    def call(env)
      @request = Rack::Request.new(env)

      auth = env['HTTP_AUTHORIZATION']

      return failure(401, 'Authorization header missing') unless auth

      token = auth[/\ABasic (.*)\z/,1]

      return failure(401, 'Authorization header malformed') unless token

      # here we simply base64-encode our user hash and pass it through
      payload = JSON.parse(Base64.decode64(token))

      env['etna.user'] = Etna::User.new(payload.map{|k,v| [k.to_sym, v]}.to_h)

      @app.call(env)
    end
  end
end
