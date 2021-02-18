require_relative '../lib/etna/client'
require_relative '../lib/etna/describe_routes'

describe Etna::Client do
  include Rack::Test::Methods

  attr_reader :app

  before(:each) do
    class Arachne
      include Etna::Application
      class Server < Etna::Server
        post '/weave/:fabric', as: :weave do
          [ 200, {}, [ @params[:fabric] ] ]
        end
        delete '/unravel/:fabric', as: :unravel do
          success_json(fabric: @params[:fabric])
        end
      end
    end
    @app = setup_app(Arachne::Server, Etna::DescribeRoutes)
    stub_request(:any, %r!^https://arachne.test/!).to_rack(@app)
  end

  after(:each) do
    Object.send(:remove_const, :Arachne)
  end

  it 'asks the etna server for routes' do
    client = Etna::Client.new('https://arachne.test', 'token')

    expect(client.routes).to eq(Arachne::Server.routes.map(&:to_hash))
  end

  it 'defines methods for handling each named route' do
    client = Etna::Client.new('https://arachne.test', 'token')
    expect(client.respond_to?(:weave)).to be_truthy

    # body form
    expect(client.weave(fabric: 'wool')).to eq('wool')

    # block form
    client.weave(fabric: 'wool') do |response|
      expect response.body to eq('wool')
    end

    # json form
    expect(client.unravel(fabric: 'wool')).to eq(fabric: 'wool')
  end

  it 'complains if the handler method is missing args' do
    client = Etna::Client.new('https://arachne.test', 'token')
    expect(client.respond_to?(:weave)).to be_truthy

    expect{client.weave}.to raise_error(ArgumentError, 'Missing required param fabric')
  end

  it 'allows users to construct an individual route' do
    client = Etna::Client.new('https://arachne.test', 'token')
    weave = client.routes.find { |e| e[:name] == 'weave' }
    expect(client.route_path(weave, {fabric: 'silk'})).to eq '/weave/silk'
  end
end
