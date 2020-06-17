describe Janus::RefreshToken do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    @user = create(
      :user,
      email: 'janus@two-faces.org',
      first_name: 'Janus', last_name: 'Bifrons'
    )
    gateway = create(:project, project_name: 'gateway', project_name_full: 'Gateway')
    @perm = create(:permission, project: gateway, user: @user, role: 'viewer')

    @token = @user.create_token!

    set_cookie([ Janus.instance.config(:token_name), @token ].join('='))
  end

  it 'updates an out-of-date token' do
    @perm.role = 'editor'
    @perm.save

    @user.refresh
    new_token = @user.create_token!

    # the token has changed
    expect(new_token).not_to eq(@token)

    # we visit janus
    auth_header(:janus)
    get('/')

    expect(last_response.status).to eq(200)

    # a new cookie is set with the new token
    cookies = parse_cookie(last_response.headers['Set-Cookie'])
    expect(cookies[Janus.instance.config(:token_name)]).to eq(new_token)
  end

  it 'ignores an up-to-date token' do
    # we visit janus
    auth_header(:janus)
    get('/')

    expect(last_response.status).to eq(200)

    # there is no cookie set
    expect(last_response.headers['Set-Cookie']).to be_nil
  end
end
