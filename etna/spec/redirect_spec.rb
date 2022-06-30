describe Etna::Redirect do
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

  def auth_header
    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'e:labors'
    ))
  end

  it "redirects to a safe url" do
    Arachne::Server.get('/test') { Etna::Redirect(@request).to('http://somewhere.example.org') }

    @app = setup_app(
      Arachne::Server,
      [ Etna::TestAuth ]
    )

    auth_header
    get('/test')

    expect(last_response.status).to eq(302)
    expect(last_response.headers['Location']).to eq(
      "https://somewhere.example.org"
    )
  end

  it "won't redirect localhost request" do
    Arachne::Server.get('/test') {
      Etna::Redirect(
        Rack::Request.new(
          @request.env.update({
            "HTTP_HOST" => "localhost",
            "SERVER_PORT" => "3000",
            "SERVER_NAME" => "localhost"
          })
        )
      ).to('http://somewhere.example.org')
    }

    @app = setup_app(
      Arachne::Server,
      [ Etna::TestAuth ]
    )

    auth_header
    get('/test')

    expect(last_response.status).to eq(422)
    expect(json_body[:errors]).to eq([ "Cannot redirect out of domain"])
  end

  it "won't redirect outside of the app domain" do
    Arachne::Server.get('/test') { Etna::Redirect(@request).to('https://somewhere.bad.org') }

    @app = setup_app(
      Arachne::Server,
      [ Etna::TestAuth ]
    )

    auth_header
    get('/test')

    expect(last_response.status).to eq(422)
    expect(json_body[:errors]).to eq([ "Cannot redirect out of domain"])
  end
end
