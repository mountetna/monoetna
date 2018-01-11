require_relative '../lib/etna/auth'
require_relative '../lib/etna/test_auth'
require 'securerandom'

describe Etna::TestAuth do
  include Rack::Test::Methods
  attr_reader :app

  before(:each) do
    class Arachne
      include Etna::Application
      class Server < Etna::Server; end
    end
  end

  after(:each) do
    Object.send(:remove_const, :Arachne)
  end

  it "creates a user without requiring validation" do
    # This exercise is the same as above
    user = nil
    Arachne::Server.get('/test') { user = @user; success('') }

    @app = setup_app(
      Arachne::Server.new(test: { }),
      [ Etna::TestAuth ]
    )
    header(*Etna::TestAuth.header(
      email: 'janus@two-faces.org',
      first: 'Janus',
      last: 'Bifrons',
      perm: 'a:labors;e:olympics,argo;v:constellations'
    ))
    get('/test')

    expect(last_response.status).to eq(200)

    expect(user).to be_a(Etna::User)
    expect(user.is_admin?('labors')).to be_truthy
    expect(user.can_edit?('constellations')).to be_falsy
  end
end

describe Etna::Auth do
  include Rack::Test::Methods
  attr_reader :app

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

  context "hmac" do
    before(:each) do
      class Arachne
        include Etna::Application
        class Server < Etna::Server; end
      end
    end

    after(:each) do
      Object.send(:remove_const, :Arachne)
    end

    it "returns an approval from an hmac" do
      # This exercise is the same as above
      Arachne::Server.get('/test') { success(nil) }

      @app = setup_app(
        Arachne::Server.new(test: {
          hmac_key: SecureRandom.hex
        }), [ Etna::Auth ]
      )

      time = DateTime.now.iso8601

      hmac = Etna::Hmac.new(
        Arachne.instance.sign,
        method: 'GET',
        host: 'example.org',
        path: '/test',
        timestamp: time,
        headers: { project_name: 'tapestry', action: 'weave' },
        id: :arachne
      )

      hmac_fields = {
        timestamp: time,
        headers: [:project_name, :action ].join(','),
        nonce: SecureRandom.hex,
        id: :arachne,
        signature: '13cf35b81cee168e50ea8521099f6dc62c4fca6c7a7b6669095b2bee71d1fcd5'
      }

      hmac_auth = hmac_fields.map{|k,v| [k,v].join('=')}.join(';')

      header('Authorization', "Hmac #{hmac_auth}")
      header('X-Etna-Project-Name', 'tapestry')
      header('X-Etna-Action', 'weave')
      header('X-Authorization-Timestamp', time)

      get('/test')
      expect(last_response.status).to eq(200)
    end
  end
end

describe Etna::Controller do
  include Rack::Test::Methods
  attr_reader :app

  before(:each) do
    class Arachne
      include Etna::Application
      class Server < Etna::Server; end
    end
  end

  after(:each) do
    Object.send(:remove_const, :Arachne)
  end

  it "applies an auth check" do
    Arachne::Server.get('/test', auth: [ :project_name, :can_edit? ] ) { success('') }

    @app = setup_app(
      Arachne::Server.new(test: { }),
      [ Etna::TestAuth ]
    )

    header(*Etna::TestAuth.header(
      email: 'janus@two-faces.org',
      perm: 'e:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(200)

    header(*Etna::TestAuth.header(
      email: 'janus@two-faces.org',
      perm: 'v:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(401)
  end

  it "applies multiple auth checks" do
    Arachne::Server.get('/test', auth: [ :project_name, [ :is_superuser?, :can_see_restricted? ] ] ) { success('') }

    @app = setup_app(
      Arachne::Server.new(test: { }),
      [ Etna::TestAuth ]
    )

    header(*Etna::TestAuth.header(
      email: 'janus@two-faces.org',
      perm: 'a:administration;V:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(200)

    header(*Etna::TestAuth.header(
      email: 'janus@two-faces.org',
      perm: 'v:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(401)
  end
end
