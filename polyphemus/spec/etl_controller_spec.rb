describe EtlController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#list' do
    it 'returns a list of etl configs for a project' do
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, secrets: { 'password' => 'shibboleth' })
      create_dummy_etl(project_name: 'athena', run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      auth_header(:editor)
      get('/api/etl/labors/configs')

      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(1)
      expect(json_body.first.keys).to match_array([
        :project_name, :etl, :name, :ran_at, :run_interval, :archived, :status, :updated_at, :created_at, :config, :comment, :secrets, :params
      ])
      expect(json_body.first[:project_name]).to eq('labors')
      expect(json_body.first[:secrets]).to eq(password: '***')
    end
  end

  context '#output' do
    it 'returns output for an etl_config' do
      output = 'A serious error happened'
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, output: output)
      auth_header(:editor)
      get(URI.encode('/api/etl/labors/output/Dummy ETL'))

      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(1)
      expect(json_body).to eq(output: output)
    end
  end

  context '#jobs' do
    it 'returns a list of etl jobs for a project' do
      auth_header(:editor)
      get('/api/etl/jobs')

      expect(last_response.status).to eq(200)
      expect(json_body.map(&:keys)).to all(match_array([:name, :schema, :secrets, :params]))
    end
  end

  context '#update' do
    before(:all) do
      create_dummy_job
    end

    after(:all) do
      remove_dummy_job
    end

    it 'updates an etl' do
      etl = create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER)

      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      expect(last_response.status).to eq(200)

      etl.refresh

      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_ONCE)
    end

    it 'clears up duplicate etls' do
      etl = create_dummy_etl(updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      etl2 = create_dummy_etl(updated_at: DateTime.now + 20, run_interval: Polyphemus::EtlConfig::RUN_NEVER)

      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      expect(last_response.status).to eq(200)

      etl.refresh
      etl2.refresh

      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.archived).to be_truthy
      expect(etl2.run_interval).to eq(Polyphemus::EtlConfig::RUN_ONCE)
      expect(etl2.archived).to be_falsy
    end

    it 'archives updated configs' do
      etl = create_dummy_etl(updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      new_config = { 'foo' => 2 }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', config: new_config)

      expect(last_response.status).to eq(200)

      etl.refresh
      expect(etl.archived).to be_truthy
      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.config).to eq('foo' => 2)

      expect(Polyphemus::EtlConfig.count).to eq(2)

      etl2 = Polyphemus::EtlConfig.last

      expect(etl2.archived).to be_falsy
      expect(etl2.config).to eq(new_config)
    end

    it 'rejects an invalid config' do
      etl = create_dummy_etl(updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      new_config = { 'blah' => 'blah' }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', config: new_config)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid configuration for etl "dummy"')
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.archived).to be_falsy
      expect(etl.config.to_h).to eq('foo' => 2)
    end

    it 'updates secrets' do
      etl = create_dummy_etl(
        updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE,
        secrets: { 'password' => 'shibboleth', 'rumor' => 'King Midas has the ears of an ass' }
      )

      new_secret = { 'rumor' => 'Midas has the ears of an ass' }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', secrets: new_secret)

      expect(last_response.status).to eq(200)
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.secrets).to eq('rumor' => 'Midas has the ears of an ass', 'password' => 'shibboleth')
    end

    it 'complains about unknown secrets' do
      etl = create_dummy_etl(updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      new_secret = { 'barber' => 'Midas has the ears of an ass' }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', secrets: new_secret)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Secrets for dummy jobs must be one of: rumor, password')
      etl.refresh
      expect(etl.secrets).to eq({})
    end

    it 'updates params' do
      etl = create_dummy_etl(
        updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE,
        params: { problem: 'present' }
      )

      new_params = { whippit: true }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', params: new_params)

      expect(last_response.status).to eq(200)
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.params.to_h).to eq({'problem' => 'present', 'whippit' => true})
    end

    it 'complains about invalid params' do
      etl = create_dummy_etl(
        updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE,
        params: {}
      )

      new_params = { problem: 'bogus', zippit: true, select_one: ["all"] }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', params: new_params)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Problem must be in: present, absent; no such param zippit; select_one must be a comma-separated string')
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.params.to_h).to eq({})
    end
  end

  context '#create' do
    before(:all) do
      create_dummy_job
    end

    after(:all) do
      remove_dummy_job
    end

    it 'creates an etl' do
      auth_header(:editor)
      json_post('/api/etl/labors/create/Dummy ETL', job_type: 'dummy')

      expect(last_response.status).to eq(200)

      expect(Polyphemus::EtlConfig.count).to eq(1)

      etl = Polyphemus::EtlConfig.last

      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.name).to eq('Dummy ETL')
    end

    it 'doesn\'t create duplicate etls' do
      etl = create_dummy_etl(updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_NEVER)

      auth_header(:editor)
      json_post('/api/etl/labors/create/Dummy ETL', job_type: 'dummy')

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('There is already an etl Dummy ETL configured for project labors')

      expect(Polyphemus::EtlConfig.count).to eq(1)
    end

    it 'rejects invalid job types' do
      auth_header(:editor)
      json_post('/api/etl/labors/create/Dummy ETL', job_type: 'mummy')

      expect(Polyphemus::EtlConfig.count).to eq(0)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('There is no such job type mummy')
    end
  end

  context '#revisions' do
    it 'returns previous revisions for an etl' do
      etl = create_dummy_etl(config: {}, updated_at: DateTime.parse('1000-01-01'), archived: true, comment: 'first pass')
      etl = create_dummy_etl(config: { 'foo': 'bar'}, updated_at: DateTime.parse('1000-02-01'), archived: true, comment: 'some tweaks')
      etl = create_dummy_etl(config: { 'foo': 'baz'}, updated_at: DateTime.parse('1000-03-01'), archived: true, comment: 'almost got it')
      etl = create_dummy_etl(config: { 'foo': 1 }, updated_at: DateTime.parse('1000-04-01'), archived: false, comment: 'final version')

      auth_header(:editor)
      get(URI.encode('/api/etl/labors/revisions/Dummy ETL'))

      expect(last_response.status).to eq(200)
      expect(json_body.map{|r| r[:config]}).to eq([
        { foo: 1 }, { foo: 'baz' }, { foo: 'bar' }, {}
      ])
      expect(json_body.map(&:keys)).to all(match_array([ :config, :updated_at, :comment ]))
    end
  end
end
