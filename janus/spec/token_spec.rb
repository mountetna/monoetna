

describe "Token Generation" do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  def auth_header(token)
    header('Authorization', "Signed-Nonce #{token}")
  end

  def request_token rsa_key, nonce, email
    payload = [
                nonce,
                Base64.strict_encode64(email)
              ].join('.')

    # sign it and base64 encode
    my_sig = Base64.strict_encode64(
      rsa_key.sign(
        OpenSSL::Digest::SHA256.new,
        payload
      )
    )

    return [payload, my_sig].join('.')
  end

  before(:each) do
    @rsa_key = OpenSSL::PKey::RSA.generate 1024

    @user = create(
      :user,
      email: 'janus@mount.etna',
      name: 'Janus',
      public_key: @rsa_key.public_key.to_s
    )
  end

  it "returns a signed timestamp token (nonce)" do
    Timecop.freeze(Time.now)
    get('/api/tokens/nonce')

    expect(last_response.status).to eq(200)
    date, sig = Base64.decode64(last_response.body).split(/\./)

    date = DateTime.parse(date)
    expect(date.iso8601).to eq(DateTime.now.iso8601)
    expect(sig).to match(/\A[0-9a-f]+\z/)
    Timecop.return
  end

  it "allows token generation with a user signature" do
    # First the user gets a nonce
    get('/api/tokens/nonce')

    # Encode your authorization request
    request = request_token(@rsa_key, last_response.body, @user.email)

    auth_header(request)
    post('/api/tokens/generate')

    expect(last_response.status).to eq(200)

    # A new token was returned
    token = last_response.body

    expect { Janus.instance.sign.jwt_decode(token) }.not_to raise_error
  end

  it "rejects token generation with an invalid nonce" do
    # This time we don't get a valid nonce from the
    # server but make up our own
    fake_nonce = Base64.strict_encode64(
      Digest::SHA256.hexdigest(
        @rsa_key.sign(OpenSSL::Digest::SHA256.new, DateTime.now.iso8601)
      )
    )

    # We sign it correctly as ourselves
    request = request_token(@rsa_key, fake_nonce, @user.email)

    header 'Authorization', "Signed-Nonce #{request}"
    post('/api/tokens/generate')

    # The request is rejected
    expect(last_response.status).to eq(401)
  end

  it "rejects token generation with an invalid email" do
    # We get the valid nonce
    get('/api/tokens/nonce')

    # But we ask for the wrong identity
    request = request_token(@rsa_key, last_response.body, 'polyphemus@mount.etna')

    header 'Authorization', "Signed-Nonce #{request}"
    post('/api/tokens/generate')

    # The request is rejected
    expect(last_response.status).to eq(401)
  end

  it "rejects token generation with an invalid signature" do
    # We get the valid nonce
    get('/api/tokens/nonce')

    # But we set a different key for signing
    bad_rsa_key = OpenSSL::PKey::RSA.generate(1024)

    request = request_token(bad_rsa_key, last_response.body, @user.email)

    header 'Authorization', "Signed-Nonce #{request}"
    post('/api/tokens/generate')

    # The request is rejected
    expect(last_response.status).to eq(401)
  end

  it "rejects token generation if the nonce is stale" do
    # First the user gets a nonce
    now = Time.now
    Timecop.freeze(now)

    get('/api/tokens/nonce')

    # Encode your authorization request
    request = request_token(@rsa_key, last_response.body, @user.email)

    # Before 60 seconds it is valid
    Timecop.freeze(now + 40) do
      auth_header(request)
      post('/api/tokens/generate')
      expect(last_response.status).to eq(200)
    end

    # After 60 seconds it is invalid
    Timecop.freeze(now + 60) do
      auth_header(request)
      post('/api/tokens/generate')
      expect(last_response.status).to eq(401)
    end

    Timecop.return
  end

  context 'task tokens' do
    before(:each) do
      @user = create(:user, name: 'Zeus Almighty', email: 'zeus@olympus.org')
      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      tannel = create(:project, project_name: 'tannel', project_name_full: 'Tunnel with an a')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')

      perm = create(:permission, project: tunnel, user: @user, role: 'administrator', privileged: true)
      perm = create(:permission, project: tannel, user: @user, role: 'administrator', privileged: true)
      perm = create(:permission, project: mirror, user: @user, role: 'viewer')
    end

    describe 'e2e' do
      it 'works' do
        # When re-recording, provide a TOKEN into the environment to run against your local development.
        if (tok = ENV['TOKEN'])
          header('Authorization', "Etna #{tok}")
        else
          auth_header(:viewer)
        end

        # When re-recording, you'll also have to enter the date you created the cassette
        #   to avoid token expiration on future tests.
        Timecop.freeze(DateTime.strptime("2021-05-06", "%Y-%m-%d"))

        VCR.use_cassette('task_token.e2e') do
          janus_client = Etna::Clients::Janus.new(host: 'https://janus.development.local', token: tok)
          token = janus_client.generate_token('task', project_name: 'ipi')
          janus_client = Etna::Clients::Janus.new(host: 'https://janus.development.local', token: token)
          janus_client.validate_task_token
        end

        Timecop.return
      end
    end

    it 'creates a task token' do
      Timecop.freeze

      header('Authorization', "Etna #{@user.create_token!}")
      post('/api/tokens/generate', project_name: 'tannel', token_type: 'task')
      expect(last_response.status).to eq(200)

      payload = header = nil

      expect {
        payload, header = Janus.instance.sign.jwt_decode(last_response.body)
      }.not_to raise_error

      # there is only one project on the token, with reduced permissions
      expect(payload["perm"]).to eq('E:tannel')

      expect(payload["task"]).to be_truthy

      expect(Time.at(payload["exp"]) - Time.now).to be_within(1).of(Janus.instance.config(:task_token_life))

      Timecop.return
    end

    it 'creates a task token with a signed nonce' do
      Timecop.freeze

      get('/api/tokens/nonce')

      request = request_token(@rsa_key, last_response.body, @user.email)
      header('Authorization', "Etna #{@user.create_token!}")
      post('/api/tokens/generate', project_name: 'tunnel', token_type: 'task')
      expect(last_response.status).to eq(200)

      payload = header = nil

      expect {
        payload, header = Janus.instance.sign.jwt_decode(last_response.body)
      }.not_to raise_error

      # there is only one project on the token, with reduced permissions
      expect(payload["perm"]).to eq('E:tunnel')

      expect(payload["task"]).to be_truthy

      expect(Time.at(payload["exp"]) - Time.now).to be_within(1).of(Janus.instance.config(:task_token_life))

      Timecop.return
    end

    it 'refuses to create a task token without project permission' do
      header('Authorization', "Etna #{@user.create_token!}")
      post('/api/tokens/generate', project_name: 'gateway', token_type: 'task')
      expect(last_response.status).to eq(401)
    end

    it 'allows supereditors create a task token without project permission' do
      admin = create(:project, project_name: 'administration', project_name_full: 'Administration')
      perm = create(:permission, project: admin, user: @user, role: 'editor')
      header('Authorization', "Etna #{@user.create_token!}")
      post('/api/tokens/generate', project_name: 'gateway', token_type: 'task')
      expect(last_response.status).to eq(200)
    end

    it 'validates a task token' do
      header('Authorization', "Etna #{@user.create_task_token!('tunnel')}")
      post('/api/tokens/validate_task')
      expect(last_response.status).to eq(200)
    end

    it 'rejects a regular token' do
      auth_header(:zeus)
      post('/api/tokens/validate_task')
      expect(last_response.status).to eq(401)
    end

    def update_payload(token, update)
      header, payload, sig = token.split('.')

      payload = payload.yield_self do |p|
        p = JSON.parse(Base64.decode64(p), symbolize_names: true)
        Base64.strict_encode64(p.merge(update).to_json)
      end

      return "#{header}.#{payload}.#{sig}"
    end

    it 'rejects a task token with admin rights' do
      token = update_payload(
        @user.create_task_token!('tunnel'),
        perm: 'A:tunnel'
      )
      header('Authorization', "Etna #{token}")
      post('/api/tokens/validate_task')
      expect(last_response.status).to eq(401)
    end

    it 'rejects a task token for the administrator group' do
      token = update_payload(
        @user.create_task_token!('tunnel'),
        perm: 'E:administration'
      )
      header('Authorization', "Etna #{token}")
      post('/api/tokens/validate_task')
      expect(last_response.status).to eq(401)
    end

    it 'rejects a task token for an unauthorized project' do
      token = update_payload(
        @user.create_task_token!('tunnel'),
        perm: 'E:gateway'
      )
      header('Authorization', "Etna #{token}")
      post('/api/tokens/validate_task')
      expect(last_response.status).to eq(401)
    end

    it 'rejects a task token with several projects' do
      token = update_payload(
        @user.create_task_token!('tunnel'),
        perm: 'E:tunnel,mirror'
      )
      header('Authorization', "Etna #{token}")
      post('/api/tokens/validate_task')
      expect(last_response.status).to eq(401)
    end

    it 'rejects a task token with elevated privileges' do
      token = update_payload(
        @user.create_task_token!('mirror'),
        perm: 'E:mirror'
      )
      header('Authorization', "Etna #{token}")
      post('/api/tokens/validate_task')
      expect(last_response.status).to eq(401)
    end
  end
end
