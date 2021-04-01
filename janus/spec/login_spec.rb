describe AuthorizationController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context 'password login' do
    before(:each) do
      clear_cookies
      @refer = JANUS_URL
      @password = 'password'
      @user = create(
        :user,
        email: 'janus@two-faces.org',
        name: "Janus",
        pass_hash: Janus.instance.sign.hash_password(@password)
      )
    end

    it 'gets a simple form' do
      get("/login?refer=#{@refer}")

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/value='#{@refer}'/)
    end

    it 'redirects with a valid cookie' do
      set_cookie([ Janus.instance.config(:token_name), @user.create_token! ].join('='))

      get("/login?refer=#{@refer}")

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq(@refer)
    end

    it 'gets a simple form with an expired cookie' do
      Timecop.freeze(Time.now - Janus.instance.config(:token_life) - 10) do
        set_cookie([ Janus.instance.config(:token_name), @user.create_token! ].join('='))
      end

      get("/login?refer=#{@refer}")

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/value='#{@refer}'/)
    end

    it 'gets a simple form with out-of-date credentials' do
      set_cookie([ Janus.instance.config(:token_name), @user.create_token! ].join('='))

      # Give Janus a project
      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      perm = create(:permission, project: gateway, user: @user, role: 'editor')

      get("/login?refer=#{@refer}")

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/value='#{@refer}'/)
    end

    it 'complains without credentials' do
      json_post( 'validate-login', {} )
      expect(last_response.status).to eq(422)
    end

    it 'validates a password' do
      form_post(
        'validate-login', 
        email: @user.email,
        password: 'bassboard',
        refer: @refer
      )
      expect(last_response.status).to eq(422)
    end

    it 'sets a cookie with the token on success' do
      form_post(
        'validate-login', 
        email: @user.email,
        password: 'password',
        refer: @refer
      )
      expect(last_response.status).to eq(302)
      cookies = parse_cookie(last_response.headers['Set-Cookie'])

      # a valid token is present under the token name
      token = cookies[Janus.instance.config(:token_name)]
      expect{Janus.instance.sign.jwt_decode(token)}.not_to raise_error

      # the cookie is restricted to the token domain
      expect(cookies["domain"]).to eq(Janus.instance.config(:token_domain))

      # the cookie is secure (https only)
      expect(cookies.has_key?("secure")).to be_truthy

      # the cookie is accessible to javascript
      expect(cookies.has_key?("HttpOnly")).to be_falsy

      # the cookie is restricted to the same site
      expect(cookies["SameSite"]).to eq('Strict')
    end

    context 'cookie expiration time' do
      it 'sets the expiration time on the cookie correctly for a new token' do
        form_post(
          'validate-login', 
          email: @user.email,
          password: 'password',
          refer: @refer
        )
        cookies = parse_cookie(last_response.headers['Set-Cookie'])
        token = cookies[Janus.instance.config(:token_name)]
        cookie_time = Time.parse(cookies['expires'])
        payload = nil
        expect {
          payload, headers = Janus.instance.sign.jwt_decode(token)
        }.not_to raise_error
        expect(cookie_time).to be_within(1).of(Time.at(payload["exp"]))
        expect(last_response.status).to eq(302)
      end

      it 'sets a fresh token on login' do
        token = nil
        Timecop.freeze(Time.now - Janus.instance.config(:token_life) - 10) do
          token = @user.create_token!
        end
        set_cookie([ Janus.instance.config(:token_name), token ].join('='))
        form_post(
          'validate-login', 
          email: @user.email,
          password: 'password',
          refer: @refer
        )
        cookies = parse_cookie(last_response.headers['Set-Cookie'])
        cookie_time = Time.parse(cookies['expires'])

        expect(cookies[Janus.instance.config(:token_name)]).not_to eq(token)
        expect(last_response.status).to eq(302)
      end
    end

    it 'redirects to refer with credentials' do
      form_post(
        'validate-login', 
        email: @user.email,
        password: @password,
        refer: @refer
      )

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq(@refer)
    end

    it 'sets a cookie with credentials' do
      refer = JANUS_URL
      form_post(
        'validate-login', 
        email: @user.email,
        password: @password,
        refer: refer
      )
      expect(rack_mock_session.cookie_jar[Janus.instance.config(:token_name)]).not_to be_empty
    end
  end

  context 'shibboleth login' do
    before(:each) do
      # a kludge for changing config
      allow(Janus.instance).to receive(:config).and_call_original
      allow(Janus.instance).to receive(:config).with(:auth_method).and_return('shibboleth')

      @refer = JANUS_URL
    end

    after(:each) do
      RSpec::Mocks.space.proxy_for(Janus.instance).reset
    end

    it 'complains if there is no email' do
      get("/login?refer=#{@refer}")

      expect(last_response.status).to eq(401)
    end

    it 'complains if there is no user' do
      email = 'janus@two-faces.org'
      header('X-Shib-Attribute', email)

      get("/login?refer=#{@refer}")

      expect(last_response.status).to eq(401)
    end

    it 'creates a token and returns a user' do
      email = 'janus@two-faces.org'
      user = create(:user, email: email, name: 'Janus')
      header('X-Shib-Attribute', email)

      get("/login?refer=#{@refer}")

      expect(last_response.status).to eq(302)
      expect(last_response.headers['Location']).to eq(@refer)

      cookies = parse_cookie(last_response.headers['Set-Cookie'])
      token = cookies[Janus.instance.config(:token_name)]
      expect{Janus.instance.sign.jwt_decode(token)}.not_to raise_error
    end
  end

  context 'task tokens' do
    before(:each) do
      @user = create(:user, name: 'Zeus Almighty', email: 'zeus@olympus.org')
      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')

      perm = create(:permission, project: tunnel, user: @user, role: 'administrator', privileged: true)
      perm = create(:permission, project: mirror, user: @user, role: 'editor')
      perm = create(:permission, project: gateway, user: @user, role: 'editor')
    end

    it 'creates a task token' do
      Timecop.freeze

      auth_header(:zeus)

      post('/api/tokens/task/generate', project_name: 'tunnel')
      expect(last_response.status).to eq(200)

      payload = header = nil

      expect{
        payload, header = Janus.instance.sign.jwt_decode(json_body[:token])
      }.not_to raise_error

      # there is only one project on the token, with reduced permissions
      expect(payload["perm"]).to eq('E:tunnel')

      expect(Time.at(payload["exp"]) - Time.now).to be_within(1).of(Janus.instance.config(:task_token_life))

      Timecop.return
    end

    it 'validates a task token' do
      auth_header(:zeus)
      post('/api/tokens/task/validate', project_name: 'tunnel')
      expect(last_response.status).to eq(200)
    end

    it 'rejects a regular token' do
      auth_header(:zeus)
      post('/api/tokens/task/validate', project_name: 'tunnel')
      expect(last_response.status).to eq(401)
    end
  end
end
