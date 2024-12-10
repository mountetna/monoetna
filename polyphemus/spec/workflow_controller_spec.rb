describe WorkflowController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context 'config initialization' do
    it 'initializes a workflow' do
      auth_header(:editor)
      json_post('/api/etl/labors/create', workflow_name: 'my workflow name', workflow_type: 'metis')

      expect(last_response.status).to eq(200)
      expect(Polyphemus::Config.count).to eq(1)
      config = Polyphemus::Config.last

      expect(config.workflow_name).to eq('my workflow name')
      expect(config.workflow_type).to eq('metis')
      # expect(config.config_id).to eq(1)  This is 100 do not know why
    end

    it 'rejects invalid workflow types' do
      auth_header(:editor)
      json_post('/api/etl/labors/create', workflow_name: 'my workflow name', workflow_type: 'not a workflow')
      expect(last_response.status).to eq(422)
    end
  end

  context 'config updates' do
    it 'updates the config itself' do
      auth_header(:editor)
      json_post('/api/etl/labors/create', workflow_name: 'my workflow name', workflow_type: 'metis')
      some_config = { 'this is a config' => 'wooo!' }
      json_post("/api/etl/labors/update/#{json_body[:config_id]}", workflow_name: 'my workflow name', workflow_type: 'metis', config: some_config)
      expect(last_response.status).to eq(200)
    end

    it 'rejects an invalid config' do
      auth_header(:editor)
      json_post('/api/etl/labors/create', workflow_name: 'my workflow name', workflow_type: 'metis')
      some_config = { 'this is a config' => 'wooo!' }
      json_post("/api/etl/labors/update/#{json_body[:config_id]}", workflow_name: 'my workflow name', workflow_type: 'metis', config: some_config)
      expect(last_response.status).to eq(422)
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

    it 'ignores archived configs' do
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, secrets: { 'password' => 'shibboleth' })
      create_dummy_etl(config: { 'foo': 'bar'}, version_number: 2, comment: 'some tweaks')
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_ARCHIVED, version_number: 3, comment: 'some tweaks')
      create_dummy_etl(config_id: 2, project_name: 'athena', run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      auth_header(:superuser)
      post('/api/etl/configs')

      expect(last_response.status).to eq(200)
      expect(json_body[:configs].length).to eq(1)
      expect(json_body[:configs].first[:project_name]).to eq('athena')
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

  context '#add_output' do
    it 'sets output for an etl_config' do
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, output: 'All is well')
      output = 'A serious error occurred.'

      auth_header(:editor)
      post('/api/etl/labors/output/1', output: output)

      expect(last_response.status).to eq(200)
      expect(json_body[:output]).to eq(output)
      expect(Polyphemus::EtlConfig.first.output).to eq(output)
    end

    it 'appends output for an etl_config' do
      create_dummy_etl(run_interval: Polyphemus::EtlConfig::RUN_NEVER, output: 'All is well')
      output = 'A serious error occurred.'

      auth_header(:editor)
      post('/api/etl/labors/output/1', output: output, append: true)

      expect(last_response.status).to eq(200)
      expect(json_body[:output]).to eq('All is well' + output)
      expect(Polyphemus::EtlConfig.first.output).to eq('All is well'+output)
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
      expect(json_body.map(&:keys)).to all(match_array([ :config, :version_number, :comment, :created_at ]))
    end
  end
end
