describe StatsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  describe '#global_stats' do
    it 'gets global stats' do
      get("/stats")

      expect(last_response.status).to eq(200)

      expect(json_body).to eq("")
    end
  end

  describe '#project_stats' do
    it 'gets project stats' do
      get("/stats/projects")

      expect(json_body).to eq("")
    end
  end
end

describe ProjectController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end
  describe '#projects' do
    it 'gets project info' do
      get("/projects")

      expect(json_body).to eq("")
    end
  end
end
