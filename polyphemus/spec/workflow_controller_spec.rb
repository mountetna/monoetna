describe WorkflowController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end


  let(:run_id) { "argo_run_id_22" }

  def create_workflow_with_configs
    auth_header(:editor)
    json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
    expect(last_response.status).to eq(200)
    some_config = { "test_key" => 'update 1' }
    json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
      workflow_name: 'my workflow name',
      workflow_type: 'test workflow',
      config: some_config
    )
    expect(last_response.status).to eq(200)
  end

  context 'config creation' do
    it 'initializes a workflow' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')

      expect(last_response.status).to eq(200)
      expect(Polyphemus::Config.count).to eq(1)
      config = Polyphemus::Config.last

      expect(config.workflow_name).to eq('my workflow name')
      expect(config.workflow_type).to eq('test workflow')
      expect(config.config_id).to eq(1)
    end

    it 'rejects invalid workflow types' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'not a workflow')
      expect(last_response.status).to eq(422)
    end
  end

  context 'config updates' do

    it 'correctly updates the config and version number' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'test_key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config: some_config
      )
      expect(last_response.status).to eq(200)
      config = Polyphemus::Config.last
      # Contains the updated config
      expect(config.config).to eq(some_config)
      # Contains the same config_id
      expect(config.config_id).to eq(json_body[:config_id])
      # Contains a new version number
      expect(config.version_number).to eq(2)

      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config: another_config
      )
      expect(last_response.status).to eq(200)
      config = Polyphemus::Config.last
      expect(config.config_id).to eq(json_body[:config_id])
      expect(config.config).to eq(another_config)
      expect(config.version_number).to eq(3)

    end

    it 'rejects an invalid config' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'not a valid key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config: some_config
      )
      expect(last_response.status).to eq(422)
    end

    it 'updates secrets and not the version number' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
      expect(last_response.status).to eq(200)
      some_secret = { 'test_secret' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        secrets: some_secret
      )
      expect(last_response.status).to eq(200)
      config = Polyphemus::Config.last
      expect(config.secrets).to eq(some_secret)
      expect(config.version_number).to eq(1)

      another_secret = { 'test_secret' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        secrets: another_secret
      )
      expect(last_response.status).to eq(200)
      config = Polyphemus::Config.last
      expect(config.config_id).to eq(json_body[:config_id])
      expect(config.secrets).to eq(another_secret)
      expect(config.version_number).to eq(1)
    end

    it 'rejects unknown secrets' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
      expect(last_response.status).to eq(200)
      some_secret = { 'not_a_secret' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        secrets: some_secret
      )
      config = Polyphemus::Config.last
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Secrets for test workflow must be one of: test_secret')
      expect(config.secrets).to eq({})
    end

  end

  context 'list configs' do
    it 'returns a list of the most recent configs for a project' do
        auth_header(:editor)
        # Create one workflow config object
        json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
        expect(last_response.status).to eq(200)
        some_config = { 'test_key' => 'update 1' }
        json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
          workflow_name: 'my workflow name',
          workflow_type: 'test workflow',
          config: some_config
        ) 
        # Create another
        json_post('/api/workflows/labors/create', workflow_name: 'my 2nd workflow name', workflow_type: 'test workflow')
        expect(last_response.status).to eq(200)
        some_config = { 'test_key' => 'update 1' }
        json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
          workflow_name: 'my workflow name',
          workflow_type: 'test workflow',
          config: some_config
        )
      auth_header(:editor)
      get('/api/workflows/labors/configs')
      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(2)
    end

  end

  context 'retrieve a config' do

    before do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'test_key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config: some_config
      )
      expect(last_response.status).to eq(200)
      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config: another_config
      )
      expect(last_response.status).to eq(200)
    end

    it 'returns the most recent config when version is not specified' do
      auth_header(:editor)
      post("/api/workflows/labors/configs/#{json_body[:config_id]}")
      expect(last_response.status).to eq(200)
      expect(json_body[:version_number]).to eq(3)
    end

    it 'returns a config by id and version' do
      auth_header(:editor)
      post("/api/workflows/labors/configs/#{json_body[:config_id]}", version: 1)
      expect(last_response.status).to eq(200)
      expect(json_body[:version_number]).to eq(1)
    end
  end

  context 'list all configs' do
    it 'returns a list of the most recent configs for every project' do
      # Create one workflow config object
      config = Polyphemus::Config.create(
        project_name: 'labors',
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )
      config = Polyphemus::Config.create(
        project_name: 'athena',
        workflow_name: 'my 2nd workflow name',
        workflow_type: 'test workflow',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )
      auth_header(:superuser)
      post('/api/workflows/configs')
      expect(last_response.status).to eq(200)
      expect(json_body[:configs].length).to eq(2)
    end

    it 'restricts the list by workflow type' do
      # Create one workflow config object
      config = Polyphemus::Config.create(
        project_name: 'labors',
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )
      config = Polyphemus::Config.create( 
        project_name: 'athena',
        workflow_name: 'my 2nd workflow name',
        workflow_type: 'metis',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )
      auth_header(:superuser)
      post('/api/workflows/configs', workflow_type: 'test workflow')
      expect(last_response.status).to eq(200)
      expect(json_body[:configs].length).to eq(1)
    end

  end

  context 'workflows' do
    it 'returns a list of workflows for a project' do
      auth_header(:editor)
      get('/api/workflows')
      expect(last_response.status).to eq(200)
      expect(json_body.map(&:keys)).to all(satisfy { |v| !v.empty? && (v - [:name, :schema, :secrets, :runtime_params]).empty? })
    end
  end

context '#revisions' do
    it 'returns previous revisions for a workflow' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'test_key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config: some_config
      )
      expect(last_response.status).to eq(200)
      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test workflow',
        config: another_config
      )
      auth_header(:editor)
      get("/api/workflows/labors/revisions/#{json_body[:config_id]}")
      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(3)
      expect(json_body[0][:version_number]).to eq(3)
      expect(json_body[1][:version_number]).to eq(2)  
      expect(json_body[2][:version_number]).to eq(1)  
    end
  end

  context 'run updates' do
    before do
      create_workflow_with_configs
    end
    it 'creates the run when it does not exist' do
      auth_header(:editor)
      the_state = {"some_state" => "wooo!"}
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        state: the_state
      )
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Run.last.state.to_h).to eq(the_state)
      expect(Polyphemus::Run.last.created_at != nil)
    end

    it 'updates the run state when it exists' do
      auth_header(:editor)
      the_state = {"some_state" => "wooo!"}
      orchestrator_metadata = { "initial" => "data" }
      new_output = 'All is well'

      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        state: the_state,
      )
      expect(last_response.status).to eq(200)
      new_state = {"some_state" => "ya_hooo!"}
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        state: new_state,
        output: new_output,
        orchestrator_metadata: orchestrator_metadata
      )
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Run.last.state.to_h).to eq(new_state)
      expect(Polyphemus::Run.last.orchestrator_metadata.to_h).to eq(orchestrator_metadata)
      expect(Polyphemus::Run.last.output).to eq(new_output)
    end

    it 'appends output to the run' do
      initial_output = 'All is well'
      additional_output = 'A serious error occurred.'

      # Create run metadata with initial output
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        output: initial_output,
      )
      expect(last_response.status).to eq(200)

      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        output: additional_output,
        append_output: true
      )
      expect(last_response.status).to eq(200)
      expect(json_body[:output]).to eq(initial_output + additional_output)
    end

  end

  context 'retrieve the workflow state' do
    before do
      create_workflow_with_configs
    end
    it 'raises and error when a workflow has not been run yet' do
      auth_header(:editor)
      get("/api/workflows/labors/run/#{run_id}")
      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('No such run argo_run_id_22')
    end

    it 'returns the state for a config when a workflow has been run' do
      auth_header(:editor)
      the_state = {:some_state => "ya_hooo!"}
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        state: the_state
      )
      expect(last_response.status).to eq(200)
      get("/api/workflows/labors/run/#{run_id}")
      expect(last_response.status).to eq(200)
      expect(json_body[:state].to_h).to eq(the_state)
    end
  end

  context 'retrieve the previous workflow state' do
    before do
      create_workflow_with_configs
      # Create a new run
      auth_header(:editor)
      the_state = {:some_state => "ya_hooo!", :another_state => "woohoo!"}
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        state: the_state
      )
    end

    it 'returns the specific state when a state is specified' do
      state_to_retrieve = [:some_state]
      post("/api/workflows/labors/run/previous/#{json_body[:config_id]}",
        version_number: 2,
        state: state_to_retrieve
      )
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(:some_state => "ya_hooo!")
    end

    it 'raises an error when a state is not found' do
      post("/api/workflows/labors/run/previous/#{json_body[:config_id]}",
        version_number: 2,
        state: [:some_state_that_doesnt_exist, :another_state_that_doesnt_exist]
      )
      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq("Requested state keys some_state_that_doesnt_exist, another_state_that_doesnt_exist not found in run.state")
    end

  end

  context 'runtime configs updates' do
    before do
      create_workflow_with_configs
    end

    it 'creates new runtime config' do
      config = {:commit  => true }
      json_post("/api/workflows/labors/runtime_configs/update/#{json_body[:config_id]}",
        config: config,
        run_interval: 60
      )
      expect(last_response.status).to eq(200)
      expect(json_body[:created_at] != nil)
      expect(json_body[:config]).to eq(config)
      expect(json_body[:run_interval]).to eq(60)
    end

    it 'updates existing runtime config' do
      config = {:commit  => true }
      json_post("/api/workflows/labors/runtime_configs/update/#{json_body[:config_id]}",
        config: config,
        run_interval: 60
      )
      expect(last_response.status).to eq(200)

      # Update metadata
      new_config = {:commit  => false}
      json_post("/api/workflows/labors/runtime_configs/update/#{json_body[:config_id]}",
        config: new_config,
        run_interval: 120
      )
      expect(last_response.status).to eq(200)
      expect(json_body[:config]).to eq(new_config)
      expect(json_body[:run_interval]).to eq(120)
    end

  end

  context 'retrieving runtime configs' do
    before do
      create_workflow_with_configs
    end

    it 'raises an error when runtime config does not exist' do
      config_id_that_doesnt_exit = 10000
      get("/api/workflows/labors/runtime_configs/#{config_id_that_doesnt_exit}")
      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('No runtime config found for config_id 10000.') 
    end

    it 'retrieves existing run metadata' do
      a_config = {:commit  => true }
      json_post("/api/workflows/labors/runtime_configs/update/#{json_body[:config_id]}",
        config: a_config,
        run_interval: 120
      )
      expect(last_response.status).to eq(200)
      
      # Retrieve and verify
      get("/api/workflows/labors/runtime_configs/#{json_body[:config_id]}")
      expect(last_response.status).to eq(200)
      expect(json_body[:config]).to eq(a_config)
    end
  end

end
