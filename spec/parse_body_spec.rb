describe Etna::ParseBody do
  include Rack::Test::Methods

  attr_reader :app

  before(:each) do
    class Arachne
      include Etna::Application
      class Server < Etna::Server
        post '/silk', action: 'web#silk'
      end
    end
    class WebController < Etna::Controller
      def silk
        success(@params[:image])
      end
    end
    @app = setup_app(Arachne::Server.new(test: {}))
  end
  after(:each) do
    Object.send(:remove_const, :Arachne)
    Object.send(:remove_const, :WebController)
  end

  it "parses json into @params" do
    post(
      '/silk',
      {
        image: 'athena-is-fat.jpg'
      }.to_json,
      {
        'CONTENT_TYPE' => 'application/json'
      }
    )

    expect(last_response.status).to eq(200)
    expect(last_response.body).to eq('athena-is-fat.jpg')
  end
end
