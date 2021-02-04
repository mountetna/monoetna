describe DataController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'rejects a non-user' do
      auth_header(:non_user)
      get("/api/#{PROJECT}/workflows/umap/steps")

      expect(last_response.status).to eq(403)
    end

    it 'gets the steps for the umap workflow' do
      auth_header(:viewer)
      get("/api/#{PROJECT}/workflows/umap/steps")

      expect(last_response.status).to eq(200)

      expect(last_response.body.length > 0).to eq(true)
    end

    it 'returns 404 for a non-existent workflow' do
      auth_header(:viewer)
      get("/api/#{PROJECT}/workflows/fancy_schmancy/steps")

      expect(last_response.status).to eq(404)

      expect(json_body[:error]).to eq('No data for workflow fancy_schmancy.')
    end
  end
end
