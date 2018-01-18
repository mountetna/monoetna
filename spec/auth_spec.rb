require_relative '../lib/etna/auth'
require 'securerandom'

describe Etna::Auth do
  include Rack::Test::Methods
  attr_reader :app

  def make_hmac(params)
    Etna::Hmac.new(
      Arachne.instance.sign, {
        method: 'GET',
        host: 'example.org',
        id: :arachne
      }.merge(params)
    )
  end

  def hmac_headers(signature, fields)
    header('Authorization', "Hmac #{signature}")
    headers = fields.delete(:headers)
    fields.merge(headers).each do |header_name, value|
      header(
        "X-Etna-#{
          header_name.to_s
            .split(/_/)
            .map(&:capitalize)
            .join('-')
        }",
        value
      )
    end
    header('X-Etna-Headers', headers.keys.join(','))
  end

  context "tokens" do
    before(:each) do
      @private_key = OpenSSL::PKey::RSA.generate 1024
      class Arachne
        include Etna::Application
        class Server < Etna::Server; end
      end
    end

    after(:each) do
      Object.send(:remove_const, :Arachne)
    end

    it "returns a user from a Janus token" do
      user = nil

      Arachne::Server.get('/test') do
        user = @user
        success(nil)
      end

      rsa_key = OpenSSL::PKey::RSA.new(@private_key)
      public_key = rsa_key.public_key.to_s

      @app = setup_app(
        Arachne::Server.new(test: { rsa_private: @private_key, rsa_public: public_key, token_algo: 'RS256' }),
        [ Etna::Auth ]
      )

      token = Arachne.instance.sign.jwt_token(
        email: 'janus@two-faces.org',
        first: 'Janus',
        last: 'Bifrons',
        perm: 'a:labors;e:olympics,argo;v:constellations'
      )

      header 'Authorization', "Basic #{token}"
      get('/test')
      expect(last_response.status).to eq(200)
      expect(user).to be_a(Etna::User)
    end

    context "invalid tokens" do
      before(:each) do
        # We generate a new public key that does not match with
        # @private_key

        @invalid_public_key = OpenSSL::PKey::RSA.generate(1024).public_key
      end

      it "fails with an invalid token" do
        Arachne::Server.get('/test') { success(nil) }

        # The app is initialized with non-matching public and
        # private keys - we are using the same app for signing
        # and validation. A more likely situation is an auth
        # application that produces signed tokens (and only
        # sets rsa_private), and client applications that
        # validate those tokens (setting rsa_public).
        @app = setup_app(
          Arachne::Server.new(test: {
            rsa_private: @private_key,
            rsa_public: @invalid_public_key, 
            token_algo: 'RS256' 
          }),
          [ Etna::Auth ]
        )

        # the token is made using rsa_private, which is the
        # signing key used to generate JWTs
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          first: 'Janus',
          last: 'Bifrons',
          perm: 'a:labors;e:olympics,argo;v:constellations'
        )

        # The test user has somehow received this invalid token
        # and tries to use it for authorization
        header 'Authorization', "Basic #{token}"
        get('/test')

        # The token is rejected because the message does not
        # match
        expect(last_response.status).to eq(401)
      end

      it "redirects with an invalid token and an auth_redirect url" do
        # This exercise is the same as above
        Arachne::Server.get('/test') { success(nil) }

        # We initialize this app, this time with a auth_redirect
        @app = setup_app(
          Arachne::Server.new(test: {
            rsa_private: @private_key,
            rsa_public: @invalid_public_key, 
            token_algo: 'RS256',
            auth_redirect: "https://janus.test"

          }),
          [ Etna::Auth ]
        )
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          first: 'Janus',
          last: 'Bifrons',
          perm: 'a:labors;e:olympics,argo;v:constellations'
        )

        header 'Authorization', "Basic #{token}"
        get('/test')

        # Instead of a 401, we are redirected to
        # the auth endpoint.
        expect(last_response.status).to eq(302)
        expect(last_response.headers['Location']).to eq(
          "https://janus.test/login?refer=http%3A%2F%2Fexample.org%2Ftest"
        )
      end
    end
  end

  context "hmac approval" do
    before(:each) do
      class Arachne
        include Etna::Application
        class Server < Etna::Server; end
      end

      # A simple get route
      Arachne::Server.get('/test') { success(@params[:project_name]) }

      # Etna::Auth will check the hmac's authenticity
      @app = setup_app( Arachne::Server.new(test: { hmac_key: SecureRandom.hex }), [ Etna::Auth ])
      @time = DateTime.now.iso8601
      @nonce = SecureRandom.hex
    end

    after(:each) do
      Object.send(:remove_const, :Arachne)
    end

    it "fails without an hmac signature" do
      # if at first we don't succeed
      get('/test')
      expect(last_response.status).to eq(401)
    end

    it "fails with a non-matching path" do
      # This signs the path /nothing
      input_hmac = make_hmac(
        timestamp: @time,
        nonce: @nonce,
        path: '/nothing',
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      # we set the http authorization header, timestamp, etc.
      hmac_headers( input_hmac.signature,
        timestamp: @time,
        id: :arachne,
        nonce: @nonce,
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      # the get fails
      get('/test')
      expect(last_response.status).to eq(401)
    end

    it "succeeds with hmac headers" do
      # This signs the path /test
      input_hmac = make_hmac(
        timestamp: @time,
        nonce: @nonce,
        path: '/test',
        headers: { project_name: 'tapestry', action: 'weave' },
      )

      # we set the http authorization header, timestamp, etc.
      hmac_headers(
        input_hmac.signature,
        timestamp: @time,
        id: :arachne,
        nonce: @nonce,
        headers: { project_name: 'tapestry', action: 'weave' },
      )

      # the get succeeds
      get('/test')
      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq('tapestry')
    end

    it "succeeds with just url params" do
      input_hmac = make_hmac(
        timestamp: @time,
        nonce: @nonce,
        path: '/test',
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      params = {
        authorization: "Hmac #{input_hmac.signature}",
        timestamp: @time,
        nonce: @nonce,
        id: 'arachne',
        headers: 'project_name,action',
        project_name: 'tapestry',
        action: 'weave',
        text: input_hmac.send(:text_to_sign)
      }

      url = params.map do |name, value|
        "X-Etna-#{
          name.to_s
            .split(/_/)
            .map(&:capitalize)
            .join('-')
        }=#{
          CGI.escape(value)
        }"
      end.join('&')

      # the get succeeds
      get("/test?#{url}")
      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq('tapestry')
    end

    it "fails with missing headers" do
      # This signs the path /test
      input_hmac = make_hmac(
        timestamp: @time,
        nonce: @nonce,
        path: '/test',
        headers: { project_name: 'tapestry', action: 'weave' },
      )

      # we set the http authorization header, timestamp, etc.
      # but we leave out the 'project_name' header
      hmac_headers( input_hmac.signature,
        timestamp: @time,
        id: :arachne,
        nonce: @nonce,
        headers: { action: 'weave' }
      )

      # the get fails
      get('/test')
      expect(last_response.status).to eq(401)
    end
  end
end
