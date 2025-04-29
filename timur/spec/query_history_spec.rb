describe QueryHistoryController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    stub_event_log
  end

  def create_query_history(params)
      create(:query_history, {
        query: 'query',
        project_name: "labors",
        comment: 'comment'
      }.merge(params))
  end

  context '#create' do
    it 'saves a query to a query history record' do
      Timecop.freeze
      query = 'encodedquerystring'
      auth_header(:viewer)

      post('/api/query_history/labors/create',
        query: query,
        comment: 'test query'
      )

      expect(last_response.status).to eq(200)

      expect(json_body[:success]).to eq(true)
      expect(QueryHistory.count).to eq(1)
      expect(QueryHistory.first.to_hash).to eq(
        comment: "test query",
        created_at: Time.now.iso8601,
        id: 1,
        project_name: "labors",
        query: query,
        user: "hercules@twelve-labors.org"
      )
      Timecop.return
    end
  end

  context '#list' do
    it 'lists queries for a user' do
      create_query_history(user: 'hera@olympus.org')
      create_query_history(user: 'hercules@twelve-labors.org')
      create_query_history(user: 'hercules@twelve-labors.org')

      auth_header(:viewer)
      get('/api/query_history/labors')

      expect(last_response.status).to eq(200)

      expect(json_body[:queries].count).to eq(2)
      expect(json_body[:queries]).to all(include(
        user: 'hercules@twelve-labors.org'
      ))
    end
  end

  context '#remove' do
    it 'removes queries for a user' do
      create_query_history(user: 'hera@olympus.org')
      create_query_history(user: 'hercules@twelve-labors.org')
      q = create_query_history(user: 'hercules@twelve-labors.org')

      auth_header(:viewer)
      delete("/api/query_history/labors/#{q.id}")

      expect(last_response.status).to eq(200)

      expect(json_body[:success]).to eq(true)
      expect(QueryHistory.count).to eq(2)
      expect(QueryHistory.select_map(:id)).not_to include(q.id)
    end

    it 'does not remove queries for a non-user' do
      create_query_history(user: 'hera@olympus.org')
      create_query_history(user: 'hercules@twelve-labors.org')
      q = create_query_history(user: 'hercules@twelve-labors.org')

      auth_header(:admin)
      delete("/api/query_history/labors/#{q.id}")

      expect(last_response.status).to eq(403)

      expect(QueryHistory.count).to eq(3)
      expect(QueryHistory.select_map(:id)).to include(q.id)
    end
  end
end
