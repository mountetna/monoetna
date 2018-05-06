describe Etna::CrossOrigin do
  include Rack::Test::Methods

  attr_reader :app
  before(:each) do
    class Arachne
      include Etna::Application
      class Server < Etna::Server
        get '/silk' do success('ok') end
      end
    end
    @app = setup_app(Arachne::Server.new(test: {}), [Etna::CrossOrigin])
  end

  after(:each) do
    Object.send(:remove_const, :Arachne)
  end

  it 'answers CORS preflight requests from the same domain' do
    subdomain = 'https://subdomain.example.org'
    header('Origin', subdomain)
    header('Access-Control-Request-Headers', 'content-type')
    options '/silk'
    expect(last_response).to be_ok
    expect(last_response.headers['Access-Control-Allow-Methods']).to eq('GET, POST, PUT, DELETE, OPTIONS')
    expect(last_response.headers['Access-Control-Allow-Origin']).to eq(subdomain)
    expect(last_response.headers['Access-Control-Allow-Headers']).to eq('content-type')
    expect(last_response.headers['Access-Control-Allow-Credentials']).to eq('true')
  end

  it 'refuses CORS preflight requests from another domain' do
    subdomain = 'https://subdomain.other.org'
    header('Origin', subdomain)
    header('Access-Control-Request-Headers', 'content-type')
    options '/silk'
    expect(last_response).to be_ok
    expect(last_response.headers['Access-Control-Allow-Methods']).to be_nil
    expect(last_response.headers['Access-Control-Allow-Origin']).to be_nil
    expect(last_response.headers['Access-Control-Allow-Headers']).to be_nil
    expect(last_response.headers['Access-Control-Allow-Credentials']).to be_nil
  end

  it 'gives an actual cross-origin response' do
    header('Origin', 'https://subdomain.example.org')
    get '/silk'
    expect(last_response.status).to eq(200)
    expect(last_response.headers['Access-Control-Allow-Origin']).to eq('https://subdomain.example.org')
    expect(last_response.headers['Access-Control-Allow-Credentials']).to eq('true')
  end
end
