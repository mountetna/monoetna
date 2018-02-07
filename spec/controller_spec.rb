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

  it "requires parameters to be set" do
    Arachne::Server.get('/test') { require_param(:project_name); success('') }

    @app = setup_app(
      Arachne::Server.new(test: { }),
      [ Etna::TestAuth ]
    )

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'e:labors'
    ))

    get('/test')
    expect(last_response.status).to eq(422)

    get('/test?project_name=labors')
    expect(last_response.status).to eq(200)
  end
end
