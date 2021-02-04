require_relative '../lib/server/controllers/workflows_controller'

describe WorkflowsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'gets a list of workflows' do
      auth_header(:viewer)
      get("/api/#{PROJECT}/workflows")

      expect(last_response.status).to eq(200)
      expect(last_response.body).to match(/umap/)
    end

    it 'rejects a non-user' do
      auth_header(:non_user)
      get("/api/#{PROJECT}/workflows")

      expect(last_response.status).to eq(403)
    end
  end
end
