describe StepsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'returns 422 for a non-user' do
      auth_header(:non_user)
      get("/api/workflows/umap/steps")

      expect(last_response.status).to eq(422)
    end

    it 'gets the steps for the umap workflow' do
      auth_header(:viewer)
      get("/api/workflows/umap/steps")

      expect(last_response.status).to eq(200)

      expect(
        json_body[:document][:tabs][0][:panes][0][:items].map{|s| s[:name].to_sym }
      ).to eq([:weight, :size, :odor])
    end

    it 'returns 404 for a non-existent workflow' do
      auth_header(:viewer)
      get("/api/workflows/fancy_schmancy/steps")

      expect(last_response.status).to eq(404)

      expect(json_body[:error]).to eq('No steps for workflow fancy_schmancy.')
    end
  end
end
