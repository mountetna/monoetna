describe EtlController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:all) do
    create_dummy_job
  end

  after(:all) do
    remove_dummy_job
  end

  context '#list' do
    it 'returns a list of etl configs for a project' do
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, secrets: { 'password' => 'shibboleth' })
      create_dummy_etl(config_id: 2, project_name: 'athena', run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      auth_header(:editor)
      get('/api/etl/labors/configs')

      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(1)
      expect(json_body.first.keys).to match_array([
        :project_name, :etl, :name, :ran_at, :run_interval, :status, :config_id, :version_number, :created_at, :config, :comment, :secrets, :params
      ])
      expect(json_body.first[:project_name]).to eq('labors')
      expect(json_body.first[:secrets]).to eq(password: '***')
    end
  end

  context '#list_all' do
    it 'returns a list of etl configs for every project' do
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, secrets: { 'password' => 'shibboleth' })
      create_dummy_etl(config_id: 2, project_name: 'athena', run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      create_dummy_etl(config_id: 3, etl: 'metis', run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      auth_header(:superuser)
      post('/api/etl/configs')

      expect(last_response.status).to eq(200)
      expect(json_body[:configs].length).to eq(3)
      expect(json_body[:configs].map(&:keys)).to all(match_array([
        :project_name, :etl, :name, :ran_at, :run_interval, :status, :config_id, :version_number, :created_at, :config, :comment, :secrets, :params
      ]))
    end

    it 'restricts the list by type' do
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, secrets: { 'password' => 'shibboleth' })
      create_dummy_etl(config_id: 2, project_name: 'athena', run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      create_dummy_etl(config_id: 3, etl: 'metis', run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      auth_header(:superuser)
      post('/api/etl/configs', job_type: 'dummy')

      expect(last_response.status).to eq(200)
      expect(json_body[:configs].length).to eq(2)
      expect(json_body[:configs].map(&:keys)).to all(match_array([
        :project_name, :etl, :name, :ran_at, :run_interval, :status, :config_id, :version_number, :created_at, :config, :comment, :secrets, :params
      ]))
    end
  end

  context '#output' do
    it 'returns output for an etl_config' do
      output = 'A serious error happened'
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, output: output)
      auth_header(:editor)
      get('/api/etl/labors/output/1')

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
      expect(json_body.map(&:keys)).to all(satisfy { |v| !v.empty? && (v - [:name, :schema, :secrets, :params]).empty? })
    end
  end

  context '#update' do
    it 'updates an etl' do
      etl = create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER)

      auth_header(:editor)
      json_post('/api/etl/labors/update/1', run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      expect(last_response.status).to eq(200)

      etl.refresh

      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_ONCE)
    end

    it 'rejects an invalid config' do
      etl = create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      new_config = { 'blah' => 'blah' }
      auth_header(:editor)
      json_post('/api/etl/labors/update/1', config: new_config)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq("Invalid configuration for etl \"dummy\"\nroot is missing required keys: foo\nproperty '/blah' is invalid: error_type=schema")
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.config.to_h).to eq('foo' => 2)
    end

    it 'updates secrets' do
      etl = create_dummy_etl(
        run_interval: Polyphemus::EtlConfig::RUN_ONCE,
        secrets: { 'password' => 'shibboleth', 'rumor' => 'King Midas has the ears of an ass' }
      )

      new_secret = { 'rumor' => 'Midas has the ears of an ass' }
      auth_header(:editor)
      json_post('/api/etl/labors/update/1', secrets: new_secret)

      expect(last_response.status).to eq(200)
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.secrets).to eq('rumor' => 'Midas has the ears of an ass', 'password' => 'shibboleth')
    end

    it 'complains about unknown secrets' do
      etl = create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      new_secret = { 'barber' => 'Midas has the ears of an ass' }
      auth_header(:editor)
      json_post('/api/etl/labors/update/1', secrets: new_secret)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Secrets for dummy jobs must be one of: rumor, password')
      etl.refresh
      expect(etl.secrets).to eq({})
    end

    it 'updates params' do
      etl = create_dummy_etl(
        run_interval: Polyphemus::EtlConfig::RUN_ONCE,
        params: { problem: 'present' }
      )

      new_params = { whippit: true }
      auth_header(:editor)
      json_post('/api/etl/labors/update/1', params: new_params)

      expect(last_response.status).to eq(200)
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.params.to_h).to eq({'problem' => 'present', 'whippit' => true})
    end

    it 'complains about invalid params' do
      etl = create_dummy_etl(
        run_interval: Polyphemus::EtlConfig::RUN_ONCE,
        params: {}
      )

      new_params = { problem: 'bogus', zippit: true, select_one: ["all"] }
      auth_header(:editor)
      json_post('/api/etl/labors/update/1', params: new_params)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Problem must be in: present, absent; no such param zippit; select_one must be a comma-separated string')
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.params.to_h).to eq({})
    end
  end

  context '#create' do
    it 'creates an etl' do
      auth_header(:editor)
      json_post('/api/etl/labors/create', name: 'Dummy ETL', job_type: 'dummy')

      expect(last_response.status).to eq(200)

      expect(Polyphemus::EtlConfig.count).to eq(1)

      etl = Polyphemus::EtlConfig.last

      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.name).to eq('Dummy ETL')
      expect(etl.config_id).to eq(1)
    end

    it 'rejects invalid job types' do
      auth_header(:editor)
      json_post('/api/etl/labors/create', name: 'Dummy ETL', job_type: 'mummy')

      expect(Polyphemus::EtlConfig.count).to eq(0)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('There is no such job type mummy')
    end
  end

  context '#revisions' do
    it 'returns previous revisions for an etl' do
      etl = create_dummy_etl(config: {}, version_number: 1, comment: 'first pass')
      etl = create_dummy_etl(config: { 'foo': 'bar'}, version_number: 2, comment: 'some tweaks')
      etl = create_dummy_etl(config: { 'foo': 'baz'}, version_number: 3, comment: 'almost got it')
      etl = create_dummy_etl(config: { 'foo': 1 }, version_number: 4, comment: 'final version')

      auth_header(:editor)
      get(URI.encode('/api/etl/labors/revisions/1'))

      expect(last_response.status).to eq(200)
      expect(json_body.map{|r| r[:config]}).to eq([
        { foo: 1 }, { foo: 'baz' }, { foo: 'bar' }, {}
      ])
      expect(json_body.map(&:keys)).to all(match_array([ :config, :version_number, :comment ]))
    end
  end
end
