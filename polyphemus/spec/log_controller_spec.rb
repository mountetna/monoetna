describe LogController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  def create_log(params)
    create(:log, {
      application: 'polyphemus',
      project_name: 'labors',
      user: 'eurystheus@twelve-labors.org',
      event: 'update',
      created_at: DateTime.now
    }.merge(params))
  end

  context '#write' do
    it 'creates a log' do
      msg = 'The Nemean Lion was slain'
      hmac_header
      json_post(
        '/api/log/labors/write',
        user: 'eurystheus@twelve-labors.org',
        event: 'update',
        message: msg,
        payload: {
          victims: 122
        }
      )
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Log.count).to eq(1)
      expect(Polyphemus::Log.first.message).to eq(msg)
    end

    it 'consolidates logs with matching entries from the past day' do
      msg = 'The Nemean Lion was slain'
      create_log(message: msg, payload: { victims: [ 'Leon' ] })
      hmac_header
      json_post(
        '/api/log/labors/write',
        user: 'eurystheus@twelve-labors.org',
        event: 'update',
        message: msg,
        payload: {
          victims: [ 'Outis', 'Aisopos' ]
        },
        consolidate: true
      )
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Log.count).to eq(1)
      expect(Polyphemus::Log.first.payload['victims']).to eq(['Leon', 'Outis', 'Aisopos'])
    end
  end

  context '#read' do
    it 'returns logs if they have restricted access' do
      create_log(message: 'The Nemean Lion was slain.')
      create_log(message: 'The Ceryneian Hind was captured.')

      auth_header(:privileged_editor)
      post('/api/log/labors/read')
      expect(last_response.status).to eq(200)
      expect(json_body[:logs].map(&:keys)).to all(contain_exactly(
        :application, :created_at, :event, :message, :project_name, :user, :id, :payload
      ))
    end

    it 'returns subsets of logs if queried' do
      create_log(message: 'The Nemean Lion was slain.')
      create_log(message: 'The Ceryneian Hind was captured.')

      auth_header(:privileged_editor)
      post('/api/log/labors/read', message: 'Ceryneian', event: '(up|down)date')
      expect(last_response.status).to eq(200)
      expect(json_body[:logs].count).to eq(1)
      expect(json_body[:logs].first[:message]).to match(/Ceryneian/)
    end
  end
  context '#payload' do
    it 'returns nothing without restricted access' do
      log = create_log(message: 'The Nemean Lion was slain.', payload: { victims: 122 })

      auth_header(:viewer)
      post("/api/log/labors/payload/#{log.id}")

      expect(last_response.status).to eq(403)
    end

    it 'returns a payload' do
      log = create_log(message: 'The Nemean Lion was slain.', payload: { victims: 122 })

      auth_header(:privileged_editor)
      post("/api/log/labors/payload/#{log.id}")

      expect(last_response.status).to eq(200)
      expect(json_body[:payload]).to eq(victims: 122)
    end

    it 'complains if the log id is not from project' do
      log = create_log(message: 'The Nemean Lion was slain.', project_name: 'lagers', payload: { victims: 122 })

      auth_header(:privileged_editor)
      post("/api/log/labors/payload/#{log.id}")

      expect(last_response.status).to eq(403)
    end
  end

  context '#hide' do
    it 'forbids normal users from hiding' do
      log = create_log(message: 'The Nemean Lion was slain.')

      auth_header(:privileged_editor)
      post("/api/log/administration/hide/#{log.id}")

      expect(last_response.status).to eq(403)
      expect(log.hidden).to be_falsy
    end

    it 'hides the log for supereditors' do
      log = create_log(message: 'The Nemean Lion was slain.')

      auth_header(:supereditor)
      post("/api/log/administration/hide/#{log.id}")

      expect(last_response.status).to eq(200)
      log.refresh
      expect(log.hidden).to be_truthy
    end
  end
end
