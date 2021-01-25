describe WorkflowsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'gets a list of workflows' do
      auth_header(:viewer)
      get("/api/workflows")

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/CONFIG/)
      expect(last_response.body).to match(/"project_name":null/)
    end

    it 'returns 422 for a non-user' do
      auth_header(:non_user)
      get("/api/workflows")

      expect(last_response.status).to eq(422)
    end
  end
end
