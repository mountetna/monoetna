require "addressable/uri"

describe Etna::Auth do
  include Rack::Test::Methods
  attr_reader :app

  def make_hmac(params)
    Etna::Hmac.new(
      Arachne.instance, {
        method: 'GET',
        host: 'example.org',
        id: :arachne
      }.merge(params)
    )
  end

  def hmac_headers(fields)
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

  def auth_header(token)
    header('Authorization', "Etna #{token}")
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

    it 'fails without an authorization header' do
      Arachne::Server.get('/test') { success(nil) }
      @app = setup_app( Arachne::Server, [ Etna::Auth ])

      # we don't use any headers
      get('/test')

      expect(last_response.status).to eq(401)
    end

    it 'succeeds without an authorization header for noauth routes' do
      Arachne::Server.get('/test', auth: { noauth: true }) { success(nil) }
      @app = setup_app(Arachne::Server, [ Etna::Auth ])

      # we don't use any headers
      get('/test')

      expect(last_response.status).to eq(200)
    end

    context 'valid tokens' do
      authenticated_user = nil
      before(:each) do
        authenticated_user = nil
        Arachne::Server.get('/test') { authenticated_user = @user; success(@user.class) }

        @public_key = OpenSSL::PKey::RSA.new(@private_key).public_key.to_s

        @app = setup_app(
          Arachne::Server,
          [ Etna::Auth ],
          test: {
            rsa_private: @private_key,
            rsa_public: @public_key,
            token_name: 'ARACHNE_TOKEN',
            token_algo: 'RS256'
          }
        )
        clear_cookies
      end

      it "returns a user from a Janus token" do
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: 'a:labors;e:olympics,argo;v:constellations'
        )

        auth_header(token)
        get('/test')
        expect(last_response.status).to eq(200)
        expect(last_response.body).to eq('Etna::User')
        expect(authenticated_user.token).to eq(token)
      end

      it 'accepts a token via cookie' do
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: 'a:labors;e:olympics,argo;v:constellations'
        )

        set_cookie("#{Arachne.instance.config(:token_name)}=#{token}")
        get('/test')
        expect(last_response.status).to eq(200)
        expect(last_response.body).to eq('Etna::User')
        expect(authenticated_user.token).to eq(token)
      end
    end

    context "community projects" do
      authenticated_user = nil
      before(:each) do
        authenticated_user = nil
        Arachne::Server.get('/:project_name', auth: { user: { can_view?: :project_name }}) {
          authenticated_user = @user; success(@user.class)
        }
        Arachne::Server.get('/generic/route') { authenticated_user = @user; success(@user.class) }

        @public_key = OpenSSL::PKey::RSA.new(@private_key).public_key.to_s

        @app = setup_app(
          Arachne::Server,
          [ Etna::Auth ],
          test: {
            rsa_private: @private_key,
            rsa_public: @public_key,
            janus: { host: 'https://janus.example.org/' },
            token_name: 'ARACHNE_TOKEN',
            token_algo: 'RS256',
            auth_redirect: 'https://janus.example.org'
          }
        )
        clear_cookies

        stub_request(:any, /janus.example.org/).
          to_return(body: {
            projects: [{
              project_name: "test_project",
              resource: true,
              requires_agreement: true
            }]
            }.to_json
          )
      end

      it "redirects a user with no existing permissions for project" do
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: 'g:constellations',
          exp: (Time.now + 6000).to_i
        )

        auth_header(token)
        get('/test_project')
        expect(last_response.status).to eq(302)
        expect(last_response.headers['Location']).to eq(
          "https://janus.example.org/test_project/cc?refer=http%3A%2F%2Fexample.org%2Ftest_project"
        )
      end

      it "permits user with existing guest access" do
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: 'a:labors;e:olympics,argo;g:test_project'
        )

        auth_header(token)
        get('/test_project')
        expect(last_response.status).to eq(200)
        expect(last_response.body).to eq('Etna::User')
        expect(authenticated_user.token).to eq(token)
      end

      it "permits access to non-project route" do
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: 'a:labors;e:olympics,argo;g:test_project'
        )

        auth_header(token)
        get('/generic/route')
        expect(last_response.status).to eq(200)
        expect(last_response.body).to eq('Etna::User')
        expect(authenticated_user.token).to eq(token)
      end
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
          Arachne::Server,
          [ Etna::Auth ],
          test: {
            rsa_private: @private_key,
            rsa_public: @invalid_public_key,
            token_algo: 'RS256'
          },
        )

        # the token is made using rsa_private, which is the
        # signing key used to generate JWTs
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: 'a:labors;e:olympics,argo;v:constellations'
        )

        # The test user has somehow received this invalid token
        # and tries to use it for authorization
        auth_header(token)
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
          Arachne::Server,
          [ Etna::Auth ],
          test: {
            rsa_private: @private_key,
            rsa_public: @invalid_public_key,
            token_algo: 'RS256',
            auth_redirect: "https://janus.example.org"
          }
        )
        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          name: 'Janus Bifrons',
          perm: 'a:labors;e:olympics,argo;v:constellations'
        )

        auth_header(token)
        get('/test')

        # Instead of a 401, we are redirected to
        # the auth endpoint.
        expect(last_response.status).to eq(302)
        expect(last_response.headers['Location']).to eq(
          "https://janus.example.org/login?refer=http%3A%2F%2Fexample.org%2Ftest"
        )
      end
    end

    context 'long-lived tokens' do
      authenticated_user = nil
      before(:each) do
        authenticated_user = nil
        Arachne::Server.get('/test') { authenticated_user = @user; success(@user.class) }

        @public_key = OpenSSL::PKey::RSA.new(@private_key).public_key.to_s

        @app = setup_app(
          Arachne::Server,
          [ Etna::Auth ],
          test: {
            rsa_private: @private_key,
            rsa_public: @public_key,
            janus: { host: 'https://janus.example.org/' },
            token_name: 'ARACHNE_TOKEN',
            token_algo: 'RS256'
          }
        )
        clear_cookies
      end

      it "calls Janus for a task token" do
        stub_request(:any, /janus.example.org/)

        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          first: 'Janus',
          last: 'Bifrons',
          task: true,
          perm: 'e:labors',
          exp: (Time.now + 6000).to_i
        )

        auth_header(token)
        get('/test')
        expect(last_response.status).to eq(200)
        expect(last_response.body).to eq('Etna::User')
        expect(authenticated_user.token).to eq(token)
        expect(WebMock).to have_requested(:post, %r!janus.example.org/api/!)
      end

      it "rejects the request if janus rejects the task token" do
        stub_request(:any, /janus.example.org/).
          to_return(status: 401)

        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          first: 'Janus',
          last: 'Bifrons',
          task: true,
          perm: 'e:labors',
          exp: (Time.now + 6000).to_i
        )

        auth_header(token)
        get('/test')
        expect(last_response.status).to eq(401)
        expect(WebMock).to have_requested(:post, %r!janus.example.org/api/!)
      end

      it "allows the request if janus is ignored" do
        stub_request(:any, /janus.example.org/).to_return(status: 401)

        token = Arachne.instance.sign.jwt_token(
          email: 'janus@two-faces.org',
          first: 'Janus',
          last: 'Bifrons',
          task: true,
          perm: 'e:labors',
          exp: (Time.now + 6000).to_i
        )

        Arachne::Server.get('/test2', auth: { ignore_janus: true }) { success(nil) }

        auth_header(token)
        get('/test2')
        expect(last_response.status).to eq(200)
        expect(WebMock).not_to have_requested(:post, %r!janus.example.org/api/!)
      end
    end

    context 'resource projects' do
      authenticated_user = nil
      before(:each) do
        authenticated_user = nil
        Arachne::Server.get('/:project_name', auth: { user: { can_view?: :project_name }}) { authenticated_user = @user; success(@user.class) }

        @public_key = OpenSSL::PKey::RSA.new(@private_key).public_key.to_s

        @app = setup_app(
          Arachne::Server,
          [ Etna::Auth ],
          test: {
            rsa_private: @private_key,
            rsa_public: @public_key,
            janus: { host: 'https://janus.example.org/' },
            token_name: 'ARACHNE_TOKEN',
            token_algo: 'RS256'
          }
        )
        clear_cookies

        stub_request(:any, /janus.example.org\/api\/user\/projects/).to_return(body: {
          projects: [{
            project_name: 'secret',
            resource: false
          }, {
            project_name: 'public',
            resource: true
          }]
        }.to_json,
          headers: {
          'Content-Type': 'application/json'
        })
      end

      context 'with task token' do
        before(:each) do
          stub_request(:any, /janus.example.org\/api\/tokens/).to_return(status: 200)
        end

        it "allows the request if user has permissions to project" do
          token = Arachne.instance.sign.jwt_token(
            email: 'janus@two-faces.org',
            first: 'Janus',
            last: 'Bifrons',
            task: true,
            perm: 'e:labors',
            exp: (Time.now + 6000).to_i
          )

          auth_header(token)
          get('/labors')
          expect(last_response.status).to eq(200)
          expect(WebMock).to have_requested(:get, %r!janus.example.org/api/user/projects!)
        end

        it "denies the request if user does not have permission to non-resource project" do
          token = Arachne.instance.sign.jwt_token(
            email: 'janus@two-faces.org',
            first: 'Janus',
            last: 'Bifrons',
            task: true,
            perm: 'e:labors',
            exp: (Time.now + 6000).to_i
          )

          auth_header(token)
          get('/secret')
          expect(last_response.status).to eq(403)
          expect(WebMock).to have_requested(:get, %r!janus.example.org/api/user/projects!).times(2)
        end

        it "allows the request if user does not have permission to resource project" do
          token = Arachne.instance.sign.jwt_token(
            email: 'janus@two-faces.org',
            first: 'Janus',
            last: 'Bifrons',
            task: true,
            perm: 'e:labors',
            exp: (Time.now + 6000).to_i
          )

          auth_header(token)
          get('/public')
          expect(last_response.status).to eq(200)
          expect(WebMock).to have_requested(:get, %r!janus.example.org/api/user/projects!)
        end
      end

      context 'for non-task token' do
        it "allows the request if user has permissions to project" do
          token = Arachne.instance.sign.jwt_token(
            email: 'janus@two-faces.org',
            first: 'Janus',
            last: 'Bifrons',
            perm: 'e:labors',
            exp: (Time.now + 6000).to_i
          )

          auth_header(token)
          get('/labors')
          expect(last_response.status).to eq(200)
          expect(WebMock).to have_requested(:get, %r!janus.example.org/api/user/projects!)
        end

        it "denies the request if user does not have permission to non-resource project" do
          token = Arachne.instance.sign.jwt_token(
            email: 'janus@two-faces.org',
            first: 'Janus',
            last: 'Bifrons',
            perm: 'e:labors',
            exp: (Time.now + 6000).to_i
          )

          auth_header(token)
          get('/secret')
          expect(last_response.status).to eq(403)
          expect(WebMock).to have_requested(:get, %r!janus.example.org/api/user/projects!).times(2)
        end

        it "allows the request if user does not have permission to resource project" do
          token = Arachne.instance.sign.jwt_token(
            email: 'janus@two-faces.org',
            first: 'Janus',
            last: 'Bifrons',
            perm: 'e:labors',
            exp: (Time.now + 6000).to_i
          )

          auth_header(token)
          get('/public')
          expect(last_response.status).to eq(200)
          expect(WebMock).to have_requested(:get, %r!janus.example.org/api/user/projects!)
        end
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

      # A download file route
      Arachne::Server.get('/download/*file_path', {auth: {hmac: true}}) { success(@params[:file_path]) }

      # we set a key for the given id
      config = {
        arachne: {
          hmac_key: SecureRandom.hex
        },
        athena: {
          hmac_key: SecureRandom.hex
        }
      }

      # Etna::Auth will check the hmac's authenticity
      @app = setup_app(
        Arachne::Server,
        [ Etna::Auth ],
        test: config
      )

      Timecop.freeze(DateTime.now)

      @time = (DateTime.now + 10).iso8601
      @nonce = SecureRandom.hex
    end

    after(:each) do
      Object.send(:remove_const, :Arachne)
      Timecop.return
    end

    it "succeeds with hmac headers" do
      # This signs the path /test
      input_hmac = make_hmac(
        expiration: @time,
        nonce: @nonce,
        path: '/test',
        headers: { project_name: 'tapestry', action: 'weave' },
      )

      # we set the hmac headers
      hmac_headers(
        signature: input_hmac.signature,
        expiration: @time,
        id: :arachne,
        nonce: @nonce,
        headers: { project_name: 'tapestry', action: 'weave' },
      )

      # the get succeeds
      get('/test')
      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq('tapestry')
    end

    it "succeeds with empty params headers" do
      # This time we use an empty hash for headers.
      # No extra params are signed.
      # Since the path is signed, any params there will be guaranteed.
      input_hmac = make_hmac(
        expiration: @time,
        nonce: @nonce,
        path: '/test',
        headers: {},
      )

      # we set the hmac headers
      hmac_headers(
        signature: input_hmac.signature,
        expiration: @time,
        id: :arachne,
        nonce: @nonce,
        headers: {}
      )

      # the get succeeds
      get('/test')
      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq('')
    end

    it "succeeds with just url params" do
      input_hmac = make_hmac(
        expiration: @time,
        nonce: @nonce,
        path: '/test',
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      uri = URI::HTTP.build(input_hmac.url_params)

      # the get succeeds
      get("#{uri.path}?#{uri.query}")

      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq('tapestry')
    end

    it "succeeds with url params and special characters" do
      raw_path = '/download/a simple+file (with a comment).txt'
      path = raw_path.split('/').map do |c|
        Addressable::URI.encode_component(c)
      end.join('/')
      input_hmac = make_hmac(
        expiration: @time,
        nonce: @nonce,
        path: path,
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      uri = URI::HTTP.build(input_hmac.url_params)

      # the get succeeds
      get("#{uri.path}?#{uri.query}")
      puts last_response
      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq(raw_path.split('/').last)
    end

    it "fails without an hmac signature" do
      # if at first we don't succeed
      get('/test')
      expect(last_response.status).to eq(401)
    end

    it "fails with a non-matching path" do
      # This signs the path /nothing
      input_hmac = make_hmac(
        expiration: @time,
        nonce: @nonce,
        path: '/nothing',
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      # we set the hmac headers
      hmac_headers(
        signature: input_hmac.signature,
        expiration: @time,
        id: :arachne,
        nonce: @nonce,
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      # the get fails
      get('/test')
      expect(last_response.status).to eq(401)
    end

    it "fails with a non-matching but known id" do
      # This signs the path /nothing as athena
      input_hmac = make_hmac(
        expiration: @time,
        nonce: @nonce,
        id: :athena,
        path: '/nothing',
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      # we set the hmac headers
      hmac_headers(
        signature: input_hmac.signature,
        expiration: @time,
        id: :arachne,
        nonce: @nonce,
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      # the get fails
      get('/test')
      expect(last_response.status).to eq(401)
    end

    it "fails if the hmac is expired" do
      Timecop.freeze(DateTime.now)

      @time = DateTime.now + 10*60

      input_hmac = make_hmac(
        expiration: @time.iso8601,
        nonce: @nonce,
        path: '/test',
        headers: { project_name: 'tapestry', action: 'weave' }
      )

      Timecop.freeze(@time + 1 )

      params = {
        authorization: "Hmac #{input_hmac.signature}",
        expiration: @time.iso8601,
        nonce: @nonce,
        id: 'arachne',
        headers: 'project_name,action',
        project_name: 'tapestry',
        action: 'weave'
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

      expect(last_response.status).to eq(401)
    end

    it "fails with missing headers" do
      # This signs the path /test
      input_hmac = make_hmac(
        expiration: @time,
        nonce: @nonce,
        path: '/test',
        headers: { project_name: 'tapestry', action: 'weave' },
      )

      # we set the hmac headers
      # but we leave out the 'project_name' header
      hmac_headers(
        signature: input_hmac.signature,
        expiration: @time,
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
