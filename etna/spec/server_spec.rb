describe Etna::Server do
  include Rack::Test::Methods

  attr_reader :app

  before(:each) do
    class Arachne
      include Etna::Application
      class Server < Etna::Server; end
    end
    class WebController < Etna::Controller
      def silk
        success('ok')
      end
    end
  end

  after(:each) do
    Object.send(:remove_const, :Arachne)
    Object.send(:remove_const, :WebController)
  end

  it 'should allow route definitions with blocks' do
    Arachne::Server.route('GET', '/silk') do
      [ 200, {}, [ 'ok' ] ]
    end
    @app = setup_app(Arachne::Server)

    get '/silk'

    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('ok')
  end

  it 'should allow route definitions with #with' do
    Arachne::Server.with(action: 'web#silk') do
      get '/silk'
    end
    @app = setup_app(Arachne::Server)

    get '/silk'

    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('ok')
  end

  it 'should allow route definitions with actions' do
    Arachne::Server.route('GET', '/silk', action: 'web#silk')
    @app = setup_app(Arachne::Server)

    get '/silk'
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('ok')
  end

  it 'supports various method verbs' do
    Arachne::Server.get('/silk') { [ 200, {}, [ 'get' ] ] }
    Arachne::Server.post('/silk') { [ 200, {}, [ 'post' ] ] }
    Arachne::Server.put('/silk') { [ 200, {}, [ 'put' ] ] }
    Arachne::Server.delete('/silk') { [ 200, {}, [ 'delete' ] ] }
    @app = setup_app(Arachne::Server)

    get '/silk'
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('get')

    post '/silk'
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('post')

    put '/silk'
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('put')

    delete '/silk'
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('delete')
  end

  it 'matches routes regardless of initial or final /' do
    Arachne::Server.get('weaver/:name') { success(@params[:name]) }
    Arachne::Server.get('/tapestry/:name') { success(@params[:name]) }
    Arachne::Server.get('textile/:name/') { success(@params[:name]) }
    @app = setup_app(Arachne::Server)

    get '/weaver/Arachne'
    expect(last_response.body).to eq('Arachne')
    get 'weaver/Arachne'
    expect(last_response.body).to eq('Arachne')
    get '/tapestry/fable-of-the-gods'
    expect(last_response.body).to eq('fable-of-the-gods')
    get 'tapestry/fable-of-the-gods/'
    expect(last_response.body).to eq('fable-of-the-gods')
    get 'textile/silk/'
    expect(last_response.body).to eq('silk')
    get '/textile/silk'
    expect(last_response.body).to eq('silk')
  end

  it 'parses route parameters' do
    Arachne::Server.get('/silk/:thread_weight/:shape') do
      [ 200, {}, [ "#{@params[:thread_weight]}-#{@params[:shape]}" ] ]
    end
    @app = setup_app(Arachne::Server)

    get '/silk/25/octagon'

    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('25-octagon')
  end

  it 'enforces route parameters over other parameters' do
    Arachne::Server.post('/silk/:thread_weight/:shape') do
      [ 200, {}, [ "#{@params[:thread_weight]}-#{@params[:shape]}" ] ]
    end
    @app = setup_app(Arachne::Server)

    post('/silk/25/octagon', { thread_weight: 35, shape: 'dodecagon' }.to_json, 'CONTENT_TYPE' => 'application/json')
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('25-octagon')

    post('/silk/25/octagon?thread_weight=35&shape=dodecagon')
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('25-octagon')
  end

  it 'allows route globbing with slashes' do
    description = 'A tapestry/weave of the Olympians.'
    Arachne::Server.get('/silk/*description') do
      [ 200, {}, [ @params[:description] ] ]
    end
    @app = setup_app(Arachne::Server)

    get URI.encode("/silk/#{description}")

    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq(description)
  end

  it 'sets route names' do
    Arachne::Server.get('/silk/:name', as: :silk_road)

    expect(Arachne::Server.routes.first.name).to eq(:silk_road)
  end

  it 'guesses route names' do
    Arachne::Server.get('/web/:silk', action: 'web#silk')

    expect(Arachne::Server.routes.first.name).to eq(:web_silk)
  end

  it 'looks up route names' do
    Arachne::Server.get('/silk/:query', as: :silk) do
      [ 200, {}, [ route_url(:silk, @params) ] ]
    end
    @app = setup_app(Arachne::Server)

    get('/silk/tree')
    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('http://example.org/silk/tree')
  end

  it "applies an auth check" do
    Arachne::Server.get('/test', auth: { user: { can_edit?: :project_name } } ) { success('') }

    @app = setup_app(Arachne::Server, [ Etna::TestAuth ])

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'e:labors',
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(200)

    header(*Etna::TestAuth.token_header(
      email: 'janus@two-faces.org',
      perm: 'v:labors'
    ))
    get('/test?project_name=labors')

    expect(last_response.status).to eq(403)
  end

  context "applying an hmac auth check" do
    before(:each) do
      Arachne::Server.get('/test', auth: { hmac: true } ) { success(@params[:project_name]) }

      @app = setup_app( Arachne::Server, [ Etna::TestAuth ])
    end

    it "allows entry if the request is hmac-authorized" do
      header(*Etna::TestAuth.hmac_header('valid'))
      get('/test')

      expect(last_response.status).to eq(200)
    end

    it "refuses entry even for valid users without an hmac" do
      header(*Etna::TestAuth.token_header(
        email: 'janus@two-faces.org',
        perm: 'e:labors'
      ))
      get('/test?project_name=labors')

      expect(last_response.status).to eq(403)
    end
  end

  it "applies multiple auth checks" do
    Arachne::Server.get('/test', auth: { user: { is_superuser?: :project_name, can_see_restricted?: :project_name } } ) { success('') }

    @app = setup_app(Arachne::Server, [ Etna::TestAuth ])

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

    expect(last_response.status).to eq(403)
  end

  context "redacting logs" do
    it 'has no effect when not configured' do
      Arachne::Server.route('POST', '/silk/:thread_weight/:shape', action: 'web#silk')

      io = StringIO.new

      @app = setup_app(Arachne::Server, nil, {test: {log_file: io}})
      Etna::Application.instance.setup_logger

      post(
        '/silk/25/octagon',
        {
          image: 'athena-is-fat.jpg'
        }.to_json,
        {
          'CONTENT_TYPE' => 'application/json'
        }
      )

      expect(io.string.include?(
        "with params {:image=>\"athena-is-fat.jpg\", :thread_weight=>\"25\", :shape=>\"octagon\"}")).to eq(true)
      expect(last_response.status).to eq(200)
    end

    it 'removes root param values' do
      Arachne::Server.route('POST', '/silk/:thread_weight/:shape', action: 'web#silk', log_redact_keys: [ :shape, :image ])

      io = StringIO.new

      @app = setup_app(Arachne::Server, nil, {test: {log_file: io}})
      Etna::Application.instance.setup_logger

      post(
        '/silk/25/octagon',
        {
          image: 'athena-is-fat.jpg'
        }.to_json,
        {
          'CONTENT_TYPE' => 'application/json'
        }
      )

      expect(io.string.include?(
        "with params {:image=>\"*\", :thread_weight=>\"25\", :shape=>\"*\"}")).to eq(true)
      expect(last_response.status).to eq(200)
    end

    it 'removes nested values' do
      Arachne::Server.route('POST', '/silk/:thread_weight/:shape', action: 'web#silk', log_redact_keys: [ :scary ])

      io = StringIO.new

      @app = setup_app(Arachne::Server, nil, {test: {log_file: io}})
      Etna::Application.instance.setup_logger

      post(
        '/silk/25/octagon',
        {
          something: {
            deep: {
              dark: {
                scary: ["abc", "123", "BOO!"]
              }
            }
          }
        }.to_json,
        {
          'CONTENT_TYPE' => 'application/json'
        }
      )

      expect(io.string.include?(
        "with params {:something=>{:deep=>{:dark=>{:scary=>\"*\"}}}, :thread_weight=>\"25\", :shape=>\"octagon\"}")).to eq(true)
      expect(last_response.status).to eq(200)
    end
  end
end
