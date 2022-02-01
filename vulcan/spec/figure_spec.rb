require_relative '../lib/server/controllers/figure_controller'

describe FigureController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#fetch' do
    it 'returns a list of figures' do
      auth_header(:viewer)
      get("/api/athena/figures")

      expect(last_response.status).to eq(200)
      expect(json_body[:figures]).to eql(["test_concurrent_workflow.cwl"])
    end
  end
end
