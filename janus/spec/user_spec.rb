describe User do
  it 'returns a JWT' do
    user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
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
    expect(payload['name']).to eq(user.name)

    # The permissions are grouped by role (a,e,v) with upper case for privileged access
    expect(payload['perm']).to eq('V:tunnel;e:gateway,mirror')
  end

  it 'expires the JWT' do
    now = Time.now
    Timecop.freeze(now)

    user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

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

  it 'complains if the email address is not all-lower-case' do
    expect {
      create(:user, name: 'Janus Bifrons', email: 'Janus@two-faces.org')
    }.to raise_error(Sequel::ValidationFailed)
  end

  it 'sets flags on the user' do
    user = create(
      :user, name: 'Janus Bifrons', email: 'janus@two-faces.org',
      flags: [ 'doors', 'portals', 'ports' ]
    )

    # the token reports the flags
    token = user.create_token!
    payload, headers = Janus.instance.sign.jwt_decode(token)

    expect(payload['flags']).to eq(user.flags.join(';'))
  end

  context 'guest-setter' do
    it 'adds guest permission if user has agreed to coc, not on project already' do
      user = create(
        :user, name: 'Portunus', email: 'portunus@two-faces.org'
      )
      gateway = create(
        :project, project_name: 'gateway', project_name_full: 'Gateway',
        resource: true, requires_agreement: true
      )
      cc = create(
        :cc_agreement, user_email: 'portunus@two-faces.org', cc_text: 'blah blah blah',
        agreed: true, project_name: 'gateway'
      )

      user.set_guest_permissions!('gateway')
      expect(Permission.count).to eq(1)

      perm = Permission.first
      expect(perm.project_name).to eq('gateway')
      expect(perm.user).to eq(user)
      expect(perm.role).to eq('guest')
    end

    it 'does nothing if user has agreed to coc, but already has viewer+ permission' do
      user = create(
        :user, name: 'Portunus', email: 'portunus@two-faces.org'
      )
      gateway = create(
        :project, project_name: 'gateway', project_name_full: 'Gateway',
        resource: true, requires_agreement: true
      )
      permission = create(
        :permission, project: gateway, user: user, role: 'viewer'
      )
      cc = create(
        :cc_agreement, user_email: 'portunus@two-faces.org', cc_text: 'blah blah blah',
        agreed: true, project_name: 'gateway'
      )

      user.set_guest_permissions!('gateway')
      expect(Permission.count).to eq(1)

      perm = Permission.first
      expect(perm.project_name).to eq('gateway')
      expect(perm.user).to eq(user)
      expect(perm.role).to eq('viewer')
    end

    it 'does nothing if user disagreed to coc, but already has viewer+ permission' do
      user = create(
        :user, name: 'Portunus', email: 'portunus@two-faces.org'
      )
      gateway = create(
        :project, project_name: 'gateway', project_name_full: 'Gateway',
        resource: true, requires_agreement: true
      )
      permission = create(
        :permission, project: gateway, user: user, role: 'viewer'
      )
      cc = create(
        :cc_agreement, user_email: 'portunus@two-faces.org', cc_text: 'blah blah blah',
        agreed: false, project_name: 'gateway'
      )

      user.set_guest_permissions!('gateway')
      expect(Permission.count).to eq(1)

      perm = Permission.first
      expect(perm.project_name).to eq('gateway')
      expect(perm.user).to eq(user)
      expect(perm.role).to eq('viewer')
    end

    it 'removes guest permission if user has guest-level permission & most recently did not agree to coc' do
      user = create(
        :user, name: 'Portunus', email: 'portunus@two-faces.org'
      )
      gateway = create(
        :project, project_name: 'gateway', project_name_full: 'Gateway',
        resource: true, requires_agreement: true
      )
      permission = create(
        :permission, project: gateway, user: user, role: 'guest'
      )

      now = Time.now
      cc1 = create(
        :cc_agreement, user_email: 'portunus@two-faces.org', cc_text: 'blah blah blah',
        agreed: true, project_name: 'gateway',
        created_at: now-10
      )
      cc2 = create(
        :cc_agreement, user_email: 'portunus@two-faces.org', cc_text: 'blah blah blah',
        agreed: false, project_name: 'gateway',
        created_at: now
      )

      user.set_guest_permissions!('gateway')
      expect(Permission.count).to eq(0)
    end

    it 'does nothing if user has guest-level permission but never agreed or disagreed to coc' do
      user = create(
        :user, name: 'Portunus', email: 'portunus@two-faces.org'
      )
      gateway = create(
        :project, project_name: 'gateway', project_name_full: 'Gateway',
        resource: true, requires_agreement: true
      )
      permission = create(
        :permission, project: gateway, user: user, role: 'guest'
      )

      now = Time.now

      user.set_guest_permissions!('gateway')
      expect(Permission.count).to eq(1)

      perm = Permission.first
      expect(perm.project_name).to eq('gateway')
      expect(perm.user).to eq(user)
      expect(perm.role).to eq('guest')
    end
  end
end

describe UserController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#update_key' do
    it 'updates a key' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(2048)

      auth_header(:janus)
      post('/api/user/update_key', pem: rsa_key.public_key.to_s)

      user.refresh
      expect(last_response.status).to eq(200)
      expect(user.public_key).to eq(rsa_key.public_key.to_s)
    end

    it 'complains if the key is not in PEM format' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      auth_header(:janus)
      post('/api/user/update_key', pem: 'I am the very model of a modern major-general')

      # Janus complains
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Improperly formatted key')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end

    it 'complains if the key is not 2048 bits' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(1024)

      auth_header(:janus)
      post('/api/user/update_key', pem: rsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Key must be 2048 bits')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end

    it 'complains if the key is not RSA' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      dsa_key = OpenSSL::PKey::DSA.generate(2048)

      auth_header(:janus)
      post('/api/user/update_key', pem: dsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Key must be RSA')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end

    it 'complains if the key is private' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(2048)

      auth_header(:janus)
      post('/api/user/update_key', pem: rsa_key.to_s)

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
      post('/api/user/update_key', pem: rsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(403)
      expect(json_body[:error]).to eq('User not found')
    end

    it 'requires authorization' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')
      rsa_key = OpenSSL::PKey::RSA.generate(2048)

      post('/api/user/update_key', pem: rsa_key.public_key.to_s)

      # Janus complains
      expect(last_response.status).to eq(401)
      expect(last_response.body).to eq('You are unauthorized')

      # The user's key is unset
      user.refresh
      expect(user.public_key).to be_nil
    end
  end

  context '#projects' do
    it 'returns the user\'s projects' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')
      door = create(:project, project_name: 'door', project_name_full: 'Door')

      # the JWT will include a string encoding these permissions
      perm = create(:permission, project: tunnel, user: user, role: 'viewer', privileged: true)
      perm = create(:permission, project: mirror, user: user, role: 'editor')
      perm = create(:permission, project: gateway, user: user, role: 'editor')

      auth_header(:janus)
      get('/api/user/projects')

      expect(last_response.status).to eq(200)
      expect(json_body[:projects]).to eq([{
        project_name: "tunnel",
        project_name_full: "Tunnel",
        role: "viewer",
        privileged: true,
        resource: false,
        requires_agreement: false,
        cc_text: "",
      }, {
        project_name: "mirror",
        project_name_full: "Mirror",
        role: "editor",
        privileged: nil,
        resource: false,
        requires_agreement: false,
        cc_text: "",
      }, {
        project_name: "gateway",
        project_name_full: "Gateway",
        role: "editor",
        privileged: nil,
        resource: false,
        requires_agreement: false,
        cc_text: "",
      }])
    end

    it 'includes resource projects' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')
      door = create(:project, project_name: 'door', project_name_full: 'Door', resource: true)

      auth_header(:janus)
      get('/api/user/projects')

      expect(last_response.status).to eq(200)
      expect(json_body[:projects]).to eq([{
        project_name: "door",
        project_name_full: "Door",
        resource: true,
        requires_agreement: false,
        cc_text: "",
      }])
    end

    it 'includes community projects' do
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org')

      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')
      door = create(:project, project_name: 'door', project_name_full: 'Door', resource: true, requires_agreement: true)

      auth_header(:janus)
      get('/api/user/projects')

      expect(last_response.status).to eq(200)
      expect(json_body[:projects]).to eq([{
        project_name: "door",
        project_name_full: "Door",
        resource: true,
        requires_agreement: true,
        cc_text: "",
      }])
    end
  end

  context '#info' do
    it 'returns the user public key fingerprint' do
      pkey = OpenSSL::PKey::RSA.new(1024)
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org', public_key: pkey.public_key)

      auth_header(:janus)

      get('/api/user/info')

      expect(last_response.status).to eq(200)
      expect(json_body[:user]).to match(
        email: "janus@two-faces.org",
        name: "Janus Bifrons",
        public_key: user.key_fingerprint,
        joined_at: user.created_at.iso8601
      )
    end
  end

  context '#fetch_all' do
    it 'returns all users for superusers' do
      zeus = create(:user, name: 'Zeus Almighty', email: 'zeus@olympus.org')
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org', flags: ['inside', 'outside'])

      gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
      tunnel = create(:project, project_name: 'tunnel', project_name_full: 'Tunnel')
      mirror = create(:project, project_name: 'mirror', project_name_full: 'Mirror')
      door = create(:project, project_name: 'door', project_name_full: 'Door')

      # the JWT will include a string encoding these permissions
      perm = create(:permission, project: tunnel, user: user, role: 'viewer', privileged: true)
      perm = create(:permission, project: mirror, user: user, role: 'editor')
      perm = create(:permission, project: gateway, user: user, role: 'editor')

      auth_header(:zeus)

      get('/api/users')

      expect(last_response.status).to eq(200)
      expect(json_body[:users].length).to eq(2)

      janus = json_body[:users].select { |u| u[:email] == 'janus@two-faces.org' }.first
      zeus = json_body[:users].select { |u| u[:email] == 'zeus@olympus.org' }.first

      expect(janus[:flags]).to match_array(['inside', 'outside'])
      expect(janus[:projects]).to match_array(['gateway', 'tunnel', 'mirror'])

      expect(zeus[:flags]).to eq(nil)
      expect(zeus[:projects]).to eq([])
    end

    it 'rejects non-superusers' do
      user = create(:user, name: 'Zeus Almighty', email: 'zeus@olympus.org')
      user = create(:user, name: 'Janus Bifrons', email: 'janus@two-faces.org', flags: ['inside', 'outside'])

      auth_header(:janus)

      get('/api/users')

      expect(last_response.status).to eq(403)
    end
  end
end
