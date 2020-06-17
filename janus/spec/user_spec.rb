describe User do
  it 'returns a JWT' do
    user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
    gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
    tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
    mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')

    # the JWT will include a string encoding these permissions
    perm = create(:permission, project: tunnel, user: user, role: 'viewer', privileged: true)
    perm = create(:permission, project: mirror, user: user, role: 'editor')
    perm = create(:permission, project: gateway, user: user, role: 'editor')

    token = user.create_token!

    expect(token).to match(%r!^[\w\-,]+\.[\w\-,]+\.[\w\-,]+$!)

    payload, headers = Janus.instance.sign.jwt_decode(token)

    expect(payload['email']).to eq(user.email)
    expect(payload['first']).to eq(user.first_name)
    expect(payload['last']).to eq(user.last_name)

    # The permissions are grouped by role (a,e,v) with upper case for privileged access
    expect(payload['perm']).to eq('V:tunnel;e:gateway,mirror')
  end

  it 'expires the JWT' do
    now = Time.now
    Timecop.freeze(now)

    user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')

    token = user.create_token!

    rsa_private = Janus.instance.sign.private_key
    rsa_public = rsa_private.public_key

    expect {
        payload, headers = JWT.decode(
        token,
        rsa_public,
        true,
        algorithm: Janus.instance.config(:token_algo)
      )
    }.not_to raise_error(JWT::ExpiredSignature)

    Timecop.freeze(now + Janus.instance.config(:token_life) + 10) do
      expect {
        payload, headers = JWT.decode(
          token,
          rsa_public,
          true,
          algorithm: Janus.instance.config(:token_algo)
        )
      }.to raise_error(JWT::ExpiredSignature)
    end

    Timecop.return
  end
end

describe UserController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#update_key' do
    it 'updates a key' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(2048)

      auth_header(:janus)
      post('/update_key', pem: rsa_key.public_key.to_s)

      user.refresh
      expect(last_response.status).to eq(200)
      expect(user.public_key).to eq(rsa_key.public_key.to_s)
    end

    it 'complains if the key is not in PEM format' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')

      auth_header(:janus)
      post('/update_key', pem: 'I am the very model of a modern major-general')

      # Janus complains
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Improperly formatted key')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end

    it 'complains if the key is not 2048 bits' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(1024)

      auth_header(:janus)
      post('/update_key', pem: rsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Key must be 2048 bits')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end

    it 'complains if the key is not RSA' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      dsa_key = OpenSSL::PKey::DSA.generate(2048)

      auth_header(:janus)
      post('/update_key', pem: dsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Key must be RSA')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end

    it 'complains if the key is private' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(2048)

      auth_header(:janus)
      post('/update_key', pem: rsa_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Cannot use a private key')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end

    it 'complains if the user is non-existent' do
      rsa_key = OpenSSL::PKey::RSA.generate(2048)

      auth_header(:janus)
      post('/update_key', pem: rsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('User not found')
    end

    it 'requires authorization' do
      user = create(:user, first_name: 'Janus', last_name: 'Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(2048)

      post('/update_key', pem: rsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(401)
      expect(last_response.body).to eq('Authorization header missing')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end
  end
end
