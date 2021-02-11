require_relative '../lib/server/controllers/workflows_controller'

describe WorkflowsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  describe 'loading' do
    ::Dir.glob(::File.join(::File.dirname(__FILE__), '..', 'lib', 'server', 'workflows', '*.{yaml,yml,cwl}')).each do |yml|
      it "#{yml} works" do
        Etna::Workflow.from_yaml_file(File.basename(yml), File.dirname(yml))
      end
    end
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
