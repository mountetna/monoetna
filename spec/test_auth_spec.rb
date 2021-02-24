require_relative '../lib/etna/test_auth'

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
      Arachne::Server,
      [ Etna::TestAuth ]
    )

    token_params = {
      email: 'janus@two-faces.org',
      name: 'Janus Bifrons',
      perm: 'a:labors;e:olympics,argo;v:constellations'
    }

    header(*Etna::TestAuth.token_header(token_params))
    get('/test')

    expect(last_response.status).to eq(200)

    expect(user).to be_a(Etna::User)
    expect(user.is_admin?('labors')).to be_truthy
    expect(user.can_edit?('constellations')).to be_falsy
    expect(user.token).to eq("something.#{Base64.strict_encode64(token_params.to_json)}")
  end

  it 'allows noauth routes through without a header' do
    # This exercise is the same as above
    user = nil
    Arachne::Server.get('/test', auth: { noauth: true }) { success('') }

    @app = setup_app(
      Arachne::Server,
      [ Etna::TestAuth ]
    )
    get('/test')

    expect(last_response.status).to eq(200)
  end
end

