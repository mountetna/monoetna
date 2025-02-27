describe StatsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context "#file_count_by_project" do
    before(:each) do
      create_file('athena', 'wisdom.txt', WISDOM)
      create_file('labors', 'helmet.jpg', HELMET)
      create_file('labors', 'helmet-shiny.jpg', SHINY_HELMET)
    end

    it "returns file count for all projects if none specified" do
      token_header(:supereditor)
      get('/api/stats/files')

      expect(last_response.status).to eq(200)

      expect(json_body).to eq({
        athena: 1,
        labors: 2,
      })
    end

    it "returns stats only for specified projects" do
      token_header(:supereditor)
      get('/api/stats/files', projects: ['athena'])

      expect(last_response.status).to eq(200)

      expect(json_body).to eq({
        athena: 1,
      })
    end

    it "only allows supereditors" do
      token_header(:editor)
      get('/api/stats/files')

      expect(last_response.status).to eq(403)
    end
  end

  context "#byte_count_by_project" do
    before(:each) do
      create_file('athena', 'wisdom.txt', WISDOM)
      create_file('labors', 'helmet.jpg', HELMET)
      create_file('labors', 'helmet-shiny.jpg', SHINY_HELMET)
    end

    it "returns file count for all projects if none specified" do
      token_header(:supereditor)
      get('/api/stats/bytes')

      expect(last_response.status).to eq(200)

      expect(json_body).to eq({
        athena: WISDOM.bytesize,
        labors: (HELMET + SHINY_HELMET).bytesize,
      })
    end

    it "returns stats only for specified projects" do
      token_header(:supereditor)
      get('/api/stats/bytes', projects: ['athena'])

      expect(last_response.status).to eq(200)

      expect(json_body).to eq({
        athena: WISDOM.bytesize,
      })
    end

    it "only allows supereditors" do
      token_header(:editor)
      get('/api/stats/bytes')

      expect(last_response.status).to eq(403)
    end
  end

end
