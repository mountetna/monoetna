describe ClientController do
  include Rack::Test::Methods
  def app
    OUTER_APP
  end
  context '#index' do
    it 'returns the map view html' do
      token_header(:viewer)
      get('/athena')

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/CONFIG/)
    end
  end
end
