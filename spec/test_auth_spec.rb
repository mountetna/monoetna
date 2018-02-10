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
      Arachne::Server.new(test: { }),
      [ Etna::TestAuth ]
    )
    header(*Etna::TestAuth.token_header(
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

