describe EtlController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context '#list' do
    it 'returns a list of etl configs for a project' do
      create(:etl_config, project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      create(:etl_config, project_name: "athena", name: "Dummy Athena ETL", config: {}, etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      auth_header(:editor)
      get('/api/etl/labors/configs')

      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(1)
      expect(json_body.first.keys).to match_array([
        :project_name, :etl, :name, :ran_at, :run_interval, :archived, :status, :output, :updated_at, :created_at, :config
      ])
      expect(json_body.first[:project_name]).to eq('labors')
    end
  end

  context '#jobs' do
    it 'returns a list of etl jobs for a project' do
      auth_header(:editor)
      get('/api/etl/jobs')

      expect(last_response.status).to eq(200)
      expect(json_body.map(&:keys)).to all(eq([:name, :schema]))
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
      etl = create(:etl_config, project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy", run_interval: Polyphemus::EtlConfig::RUN_NEVER)

      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      expect(last_response.status).to eq(200)

      etl.refresh

      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_ONCE)
    end

    it 'clears up duplicate etls' do
      etl = create(:etl_config, project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy", updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_NEVER)
      etl2 = create(:etl_config, project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy", updated_at: DateTime.now + 20, run_interval: Polyphemus::EtlConfig::RUN_NEVER)

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
      etl = create(:etl_config, project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy", updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      new_config = { 'foo' => 2 }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', config: new_config)

      expect(last_response.status).to eq(200)

      etl.refresh
      expect(etl.archived).to be_truthy
      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.config).to eq({})

      expect(Polyphemus::EtlConfig.count).to eq(2)

      etl2 = Polyphemus::EtlConfig.last

      expect(etl2.archived).to be_falsy
      expect(etl2.config).to eq(new_config)
    end

    it 'rejects an invalid config' do
      etl = create(:etl_config, project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy", updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_ONCE)

      new_config = { 'blah' => 'blah' }
      auth_header(:editor)
      json_post('/api/etl/labors/update/Dummy ETL', config: new_config)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Invalid configuration for etl "dummy"')
      expect(Polyphemus::EtlConfig.count).to eq(1)
      etl.refresh
      expect(etl.archived).to be_falsy
      expect(etl.config).to eq({})
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
      json_post('/api/etl/labors/create/Dummy ETL', job_type: "dummy")

      expect(last_response.status).to eq(200)

      expect(Polyphemus::EtlConfig.count).to eq(1)

      etl = Polyphemus::EtlConfig.last

      expect(etl.run_interval).to eq(Polyphemus::EtlConfig::RUN_NEVER)
      expect(etl.name).to eq("Dummy ETL")
    end

    it 'doesn\'t create duplicate etls' do
      etl = create(:etl_config, project_name: "labors", name: "Dummy ETL", config: {}, etl: "dummy", updated_at: DateTime.now, run_interval: Polyphemus::EtlConfig::RUN_NEVER)

      auth_header(:editor)
      json_post('/api/etl/labors/create/Dummy ETL', job_type: "dummy")

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('There is already an etl Dummy ETL configured for project labors')

      expect(Polyphemus::EtlConfig.count).to eq(1)
    end

    it 'rejects invalid job types' do
      auth_header(:editor)
      json_post('/api/etl/labors/create/Dummy ETL', job_type: "mummy")

      expect(Polyphemus::EtlConfig.count).to eq(0)

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('There is no such job type mummy')
    end
  end
end
