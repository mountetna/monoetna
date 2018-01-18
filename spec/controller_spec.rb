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

  it "applies an auth check" do
    Arachne::Server.get('/test', auth: { user: { can_edit?: :project_name } } ) { success('') }

    @app = setup_app(
      Arachne::Server.new(test: { }),
      [ Etna::TestAuth ]
    )

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'e:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(200)

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'v:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(401)
  end

  context "applying an hmac auth check" do
    before(:each) do
      Arachne::Server.get('/test', auth: { hmac: true } ) { success(@params[:project_name]) }

      @app = setup_app(
        Arachne::Server.new(test: { }),
        [ Etna::TestAuth ]
      )
    end

    it "allows entry if the request is hmac-authorized" do
      header(*Etna::TestAuth.hmac_header(
        project_name: 'labors'
      ))
      get('/test')

      expect(last_response.status).to eq(200)
      expect(last_response.body).to eq('labors')
    end

    it "refuses entry even for valid users without an hmac" do
      header(*Etna::TestAuth.token_header(
        email: 'janus@two-faces.org',
        perm: 'e:labors'
      ))
      get('/test?project_name=labors')

      expect(last_response.status).to eq(401)
    end
  end

  it "applies multiple auth checks" do
    Arachne::Server.get('/test', auth: { user: { is_superuser?: :project_name, can_see_restricted?: :project_name } } ) { success('') }

    @app = setup_app(
      Arachne::Server.new(test: { }),
      [ Etna::TestAuth ]
    )

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'a:administration;V:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(200)

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'v:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(401)
  end
end
