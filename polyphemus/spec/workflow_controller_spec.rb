describe WorkflowController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    stub_event_log
  end

  let(:run_id) { "argo_run_id_22" }

  def create_workflow_with_configs(params={})
    create(:config, {
      project_name: 'labors',
      secrets: {},
      config_id: Polyphemus::Config.next_id,
      version_number: 2,
      workflow_name: 'my workflow name',
      workflow_type: 'test',
      config: { "test_key" => 'update 1' }
    }.merge(params))
  end

  context 'config creation' do
    it 'initializes a workflow' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test')

      expect(last_response.status).to eq(200)
      expect(Polyphemus::Config.count).to eq(1)
      config = Polyphemus::Config.last

      expect(config.workflow_name).to eq('my workflow name')
      expect(config.workflow_type).to eq('test')
      expect(config.config_id).to eq(1)
      expect(config.created_at).to be_within(3.second).of(Time.now)
      expect(config.updated_at).to be_within(3.second).of(Time.now)
    end

    it 'rejects invalid workflow types' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'not a workflow')
      expect(last_response.status).to eq(422)
    end

    it 'creates a runtime config when a workflow is created' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test')
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Config.count).to eq(1)
      expect(Polyphemus::RuntimeConfig.count).to eq(1)
      expect(Polyphemus::RuntimeConfig.last.config_id).to eq(Polyphemus::Config.last.config_id)
    end
  end

  context 'config updates' do
    let(:config_id) { create_workflow_with_configs(version_number: 1).config_id }

    it 'correctly updates the config and version number' do
      some_config = { 'test_key' => 'update 1' }

      auth_header(:editor)
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        config: some_config
      )

      expect(last_response.status).to eq(200)

      expect(Polyphemus::Config.count).to eq(2)

      first_config = Polyphemus::Config.last
      expect(first_config.config).to eq(some_config)
      expect(first_config.config_id).to eq(config_id)
      expect(first_config.version_number).to eq(2)


      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        config: another_config
      )
      expect(last_response.status).to eq(200)
      second_config = Polyphemus::Config.last
      expect(second_config.config_id).to eq(config_id)
      expect(second_config.config).to eq(another_config)
      expect(second_config.version_number).to eq(3)
      expect(second_config.updated_at).to_not eq(first_config.updated_at)
      expect(second_config.created_at).to_not eq(first_config.created_at)
    end

    it 'rejects an invalid config' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test')
      expect(last_response.status).to eq(200)

      some_config = { 'not a valid key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        config: some_config
      )
      expect(last_response.status).to eq(422)
    end

    it 'updates secrets and not the version number' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test')
      expect(last_response.status).to eq(200)
      some_secret = { 'test_secret' => 'update 1' }
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        secrets: some_secret
      )
      expect(last_response.status).to eq(200)
      the_config = Polyphemus::Config.last
      expect(the_config.secrets).to eq(some_secret)
      expect(the_config.version_number).to eq(1)

      another_secret = { 'test_secret' => 'update 2' }
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        secrets: another_secret
      )
      expect(last_response.status).to eq(200)
      updated_config = Polyphemus::Config.last
      expect(updated_config.config_id).to eq(config_id)
      expect(updated_config.secrets).to eq(another_secret)
      expect(updated_config.version_number).to eq(1)
      expect(updated_config.created_at).to eq(the_config.created_at)
      expect(updated_config.updated_at).to_not eq(the_config.updated_at)
    end

    it 'rejects unknown secrets' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test')
      expect(last_response.status).to eq(200)
      some_secret = { 'not_a_secret' => 'update 1' }
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        secrets: some_secret
      )
      config = Polyphemus::Config.last
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Secrets for test must be one of: test_secret')
      expect(config.secrets).to eq({})
    end

  end

  context 'retrieve a config' do
    let(:config_id) { create_workflow_with_configs(version_number: 1).config_id }

    before do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test')
      expect(last_response.status).to eq(200)
      some_config = { 'test_key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        config: some_config
      )
      expect(last_response.status).to eq(200)
      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{config_id}",
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        config: another_config
      )
      expect(last_response.status).to eq(200)
    end

    it 'returns the most recent config when version is not specified' do
      auth_header(:editor)
      get("/api/workflows/labors/configs/#{config_id}")
      expect(last_response.status).to eq(200)
      expect(json_body[:version_number]).to eq(3)
    end

    it 'returns a config by id and version' do
      auth_header(:editor)
      get("/api/workflows/labors/configs/#{config_id}?version=1")
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
        workflow_type: 'test',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )
      config = Polyphemus::Config.create(
        project_name: 'athena',
        workflow_name: 'my 2nd workflow name',
        workflow_type: 'test',
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
        workflow_type: 'test',
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
      post('/api/workflows/configs', workflow_type: 'test')
      expect(last_response.status).to eq(200)
      expect(json_body[:configs].length).to eq(1)
    end

  end

  context 'workflows' do
    it 'returns a list of workflows for a project' do
      auth_header(:editor)
      get('/api/workflows')
      expect(last_response.status).to eq(200)
      expect(json_body.map(&:keys)).to all(satisfy { |v| !v.empty? && (v - [:name, :schema, :secrets, :runtime_params, :workflow_path]).empty? })
    end
  end

  context '#revisions' do
    it 'returns previous revisions for a workflow' do
      config = create_config(
        workflow_name: 'my workflow name', workflow_type: 'test',
        version_number: 1
      )
      create_config(
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        version_number: 2,
        config_id: config.config_id,
        config: { 'test_key' => 'update 1' }
      )
      create_config(
        workflow_name: 'my workflow name',
        workflow_type: 'test',
        config_id: config.config_id,
        version_number: 3,
        config: { 'test_key' => 'update 2' }
      )

      auth_header(:editor)
      get("/api/workflows/labors/revisions/#{config.config_id}")

      expect(last_response.status).to eq(200)
      expect(json_body.length).to eq(3)
      expect(json_body[0][:version_number]).to eq(3)
      expect(json_body[1][:version_number]).to eq(2)  
      expect(json_body[2][:version_number]).to eq(1)  
    end
  end

  context 'run updates' do
    let(:config_id) { create_workflow_with_configs.config_id }

    it 'creates the run when it does not exist' do
      auth_header(:editor)
      the_state = {"some_state" => "wooo!"}
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: config_id,
        version_number: 2,
        name: "test-workflow-1999",
        state: the_state
      )
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Run.last.state.to_h).to eq(the_state)
      expect(Polyphemus::Run.last.created_at).to_not be_nil
      expect(Polyphemus::Run.last.name).to eq("test-workflow-1999")
    end

    it 'updates the run state when it exists' do
      auth_header(:editor)
      the_state = {"some_state" => "wooo!"}
      orchestrator_metadata = { "initial" => "data" }
      new_output = 'All is well'

      # Create initial run
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: config_id,
        version_number: 2,
        name: "test-workflow-1999",
        state: the_state
      )
      expect(last_response.status).to eq(200)
      initial_updated_at = Polyphemus::Run.last.updated_at

      # Update the run
      new_state = {"some_state" => "ya_hooo!"}
      json_post("/api/workflows/labors/run/update/#{run_id}",
        state: new_state,
        output: new_output,
        orchestrator_metadata: orchestrator_metadata
      )
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Run.last.state.to_h).to eq(new_state)
      expect(Polyphemus::Run.last.orchestrator_metadata.to_h).to eq(orchestrator_metadata)
      expect(Polyphemus::Run.last.output).to eq(new_output)
      expect(Polyphemus::Run.last.updated_at).to be > initial_updated_at
    end

    it 'appends output to the run' do
      initial_output = 'All is well'
      additional_output = 'A serious error occurred.'

      # Create run metadata with initial output
      auth_header(:editor)
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: config_id,
        version_number: 2,
        name: "test-workflow-1999",
        output: initial_output
      )
      expect(last_response.status).to eq(200)

      # Update with additional output
      auth_header(:editor)
      json_post("/api/workflows/labors/run/update/#{run_id}",
        output: additional_output,
        append_output: true
      )
      expect(last_response.status).to eq(200)
      expect(json_body[:output]).to eq(initial_output + additional_output)
    end
  end

  context 'retrieve the workflow state' do
    let(:config_id) { create_workflow_with_configs.config_id }

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
        config_id: config_id,
        version_number: 2,
        name: "test-workflow-1999",
        state: the_state
      )
      expect(last_response.status).to eq(200)
      get("/api/workflows/labors/run/#{run_id}")
      expect(last_response.status).to eq(200)
      expect(json_body[:state].to_h).to eq(the_state)
    end
  end

  context 'retrieve the previous workflow state' do
    let(:config_id) {
      config = create_workflow_with_configs(version_number: 1)
      config.config_id
    }

    def create_run(params)
      create(:run,  {
        run_id: run_id,
        name: "test-workflow-1999",
        config_id: config_id,
        version_number: 1,
        state: {},
        created_at: Time.now
      }.merge(params))
    end

    it 'returns nothing when there are no runs' do
      auth_header(:editor)
      post("/api/workflows/labors/run/previous/#{config_id}", state: [:current_labor])
      expect(last_response.status).to eq(200)
      expect(json_body).to eq({})
    end

    it 'returns the specific state when a state is specified' do
      create_run(state: { current_labor: "The Nemean Lion"})

      auth_header(:editor)
      post("/api/workflows/labors/run/previous/#{config_id}", state: [:current_labor])
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(:current_labor => "The Nemean Lion")
    end

    it 'returns all previous states available when requested' do
      create_run(state: { current_labor: "The Nemean Lion" }, created_at: Time.now - 20)
      create_run(state: { current_labor: "The Lernaean Hydra" }, created_at: Time.now, run_id: "argo_run_state_23")
      auth_header(:editor)
      post("/api/workflows/labors/run/previous/#{config_id}",
        state: [:current_labor],
        collect: true
      )
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(current_labor: [ "The Lernaean Hydra", "The Nemean Lion"])
    end

    it 'returns a previous state across version changes' do
      create_run(state: { current_labor: "The Nemean Lion"})
      create_workflow_with_configs(version_number: 2, config_id: config_id)

      auth_header(:editor)
      post("/api/workflows/labors/run/previous/#{config_id}", state: [:current_labor])
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(:current_labor => "The Nemean Lion")
    end

    it 'returns a previous state even if the last state was nil' do
      create_run(state: { current_labor: "The Nemean Lion" }, created_at: Time.now - 20)
      create_run(state: nil, created_at: Time.now, run_id: "argo_run_state_23")

      auth_header(:editor)
      post("/api/workflows/labors/run/previous/#{config_id}", state: [:current_labor])
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(:current_labor => "The Nemean Lion")
    end

    it 'returns a previous state if the last state was present but lacking the requested keys' do
      create_run(state: { current_labor: "The Nemean Lion" }, created_at: Time.now - 20)
      create_run(state: {}, created_at: Time.now, run_id: "argo_run_state_23")

      auth_header(:editor)
      post("/api/workflows/labors/run/previous/#{config_id}", state: [:current_labor])
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(:current_labor => "The Nemean Lion")
    end

    it 'returns nil when a state key is not found' do
      create_run(state: { current_labor: "The Nemean Lion" }, created_at: Time.now - 20)
      create_run(state: { current_labor: "The Lernaean Hydra" }, created_at: Time.now, run_id: "argo_run_state_23")

      auth_header(:editor)
      post("/api/workflows/labors/run/previous/#{config_id}",
        state: [:monster, :current_labor]
      )
      expect(last_response.status).to eq(200)
      expect(json_body).to eq(current_labor: "The Lernaean Hydra", monster: nil)
    end
  end

  context 'runtime configs updates' do
    let(:config_id) {
      config, *_ = create_workflow(
        workflow_name: 'my-cat-ingestion',
        run_interval: 0,
        runtime_config: {commit: true},
        run_start: '2025-01-16T12:00:00Z',
        run_stop: '2025-01-16T12:20:00Z'
      )

      config.config_id
    }

    it 'updates existing runtime config' do
      # Update metadata
      new_config = {:commit  => false}
      auth_header(:editor)
      json_post("/api/workflows/labors/runtime_configs/update/#{config_id}",
        config: new_config,
        run_interval: 120,
        disabled: true
      )
      expect(last_response.status).to eq(200)
      expect(json_body[:config]).to eq(new_config)
      expect(json_body[:run_interval]).to eq(120)
      expect(json_body[:disabled]).to eq(true)
    end

    it 'raises an error when the runtime config is invalid' do
      config = {:not_a_valid_param  => "not a boolean" }
      auth_header(:editor)
      json_post("/api/workflows/labors/runtime_configs/update/#{config_id}",
        config: config,
      )
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq("no such param not_a_valid_param")
    end

  end

  context 'retrieving runtime configs' do
    let(:config_id) { 
      config, *_ = create_workflow(
        workflow_name: 'my-cat-ingestion',
        run_interval: 0,
        runtime_config: {commit: true},
        run_start: '2025-01-16T12:00:00Z',
        run_stop: '2025-01-16T12:20:00Z'
      )

      config.config_id
    }

    it 'raises an error when runtime config does not exist' do
      config_id_that_doesnt_exit = 10000
      auth_header(:editor)
      get("/api/workflows/labors/runtime_configs/#{config_id_that_doesnt_exit}")
      expect(last_response.status).to eq(404)
      expect(json_body[:error]).to eq('No such config 10000') 
    end

    it 'retrieves existing runtime config metadata' do
      auth_header(:editor)
      get("/api/workflows/labors/runtime_configs/#{config_id}")
      expect(last_response.status).to eq(200)
      expect(json_body[:config]).to eq(commit: true)
    end
  end

  context 'running a workflow once' do
    let(:config_id) { create_workflow_with_configs.config_id }

    it 'successfully runs a workflow' do
      config = Polyphemus::Config.current.where(project_name: 'labors', config_id: config_id).first
      
      # Stub argo client command
      workflow_output = "Workflow Name: simple-two-jobs-ktm4g\nWorkflow UID: de6df0ed-b32c-4707-814f-11c323b0687b\n"
      allow(Open3).to receive(:capture3).and_return([workflow_output, "", double(success?: true)])
      auth_header(:editor)
      json_post("/api/workflows/labors/runtime_configs/run_once/#{config.config_id}", 
        workflow_type: "test",
        config: { commit: true }
      )
      expect(last_response.status).to eq(200)
      expect(Polyphemus::Run.last.name).to eq("simple-two-jobs-ktm4g")
      expect(Polyphemus::Run.last.run_id).to eq("de6df0ed-b32c-4707-814f-11c323b0687b")

      cmd = [
        "argo submit",
        TestManifest.as_json[:workflow_path].shellescape,
        "-p config_id=#{config.config_id}",
        "-p version_number=#{config.version_number}", 
        "-n argo",
        "-o yaml",
        "| grep -m2 -E 'name:|uid:'",
        "| sed 's/name:/Workflow Name:/'",
        "| sed 's/uid:/Workflow UID:/'"
      ]
      expect(Open3).to have_received(:capture3).with(cmd.join(" "))
    end

    it 'fails when workflow submission fails' do
      # If we don't mock anything this will fail since there is no workflow.yaml for test-workflow.yaml
      config = Polyphemus::Config.current.where(project_name: 'labors', config_id: config_id).first
      allow(Open3).to receive(:capture3).and_return(["", "Command failed", double(success?: false)])
      auth_header(:editor)
      json_post("/api/workflows/labors/runtime_configs/run_once/#{config.config_id}", 
        workflow_type: "test",
        config: { commit: true }
      )

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq("Failed to submit Argo workflow: Command failed")
    end

  end

  context 'status' do

    it 'returns the status of a workflow that has completed' do

      config, *_ = create_workflow(
        workflow_name: 'my-cat-ingestion',
        run_interval: 0,
        runtime_config:{commit: true},
        run_start: '2025-01-16T12:00:00Z',
        run_stop: '2025-01-16T12:20:00Z'
      )

      create_run(
        config, 2000,
        status: 'Succeeded',
        run_start: '2025-01-18T12:00:00Z',
        run_stop: '2025-01-18T12:20:00Z'
      )

      auth_header(:editor)
      get("/api/workflows/labors/status")
      expect(last_response.status).to eq(200)

      # Make sure the second version of the config is used (we use the latest run)
      expect(json_body.count).to eq(1)
      expect(json_body[0][:workflow_name]).to eq("my-cat-ingestion")
      expect(json_body[0][:workflow_type]).to eq("test")
      expect(json_body[0][:pipeline_state]).to eq("succeeded")
      # This should be the most recent run from version 2
      expect(json_body[0][:pipeline_finished_at]).to eq("2025-01-18T12:20:00Z")
    end

    it 'fetches the status of a workflow that has not completed' do

      config_2 = Polyphemus::Config.create(
        project_name: 'labors',
        workflow_name: 'my-cat-ingestion-2',
        workflow_type: 'cat-ingestion',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )

      runtime_config_2 = Polyphemus::RuntimeConfig.create(
        config_id: config_2.config_id,
        run_interval: 1000
      )

      # Create a run for the second config, but workflow is not finished (no orchestrator metadata)
      Polyphemus::Run.create(
        run_id: "hello-there-2",
        config_id: config_2.config_id,
        version_number: config_2.version_number,
        name: 'my-cat-ingestion21000',
      )
      
      # Mock the get_workflow_status method to return a running status
      workflow_json = '{"metadata":{"name":"my-cat-ingestion-2"},"status":{"phase":"Running","startedAt":"2025-01-13T20:58:09Z","finishedAt":"2025-01-15T20:58:09Z"}}'
      allow(Open3).to receive(:capture3).and_return([workflow_json, "", double(success?: true)])

      auth_header(:editor)
      get("/api/workflows/labors/status")
      expect(last_response.status).to eq(200)

      expect(json_body.count).to eq(1)
      expect(json_body[0][:workflow_name]).to eq("my-cat-ingestion-2")
      expect(json_body[0][:workflow_type]).to eq("cat-ingestion")
      expect(json_body[0][:pipeline_state]).to eq("running")
      # When argo is still running, we use the previous run's finished_at
      # In this case there is only one run, so it should be nil
      expect(json_body[0][:pipeline_finished_at]).to be_nil
    end

    it 'returns a nil status when the workflow hasnt been run yet and is not scheduled' do
      create_workflow(
        workflow_name: 'parked-workflow',
        run_interval: nil
      )

      auth_header(:editor)
      get("/api/workflows/labors/status")
      expect(last_response.status).to eq(200)
      expect(json_body.count).to eq(1)
      expect(json_body[0][:workflow_name]).to eq("parked-workflow")
      expect(json_body[0][:workflow_type]).to eq("test")
      expect(json_body[0][:pipeline_state]).to be_nil
      expect(json_body[0][:pipeline_finished_at]).to be_nil
    end

    it 'returns a pending status when the workflow hasnt been run yet and is scheduled' do
      create_workflow(
        workflow_name: 'pending-workflow',
        run_interval: 1000
      )

      auth_header(:editor)
      get("/api/workflows/labors/status")
      expect(last_response.status).to eq(200)
      expect(json_body.count).to eq(1)
      expect(json_body[0][:workflow_name]).to eq("pending-workflow")
      expect(json_body[0][:workflow_type]).to eq("test")
      expect(json_body[0][:pipeline_state]).to eq("pending")
      expect(json_body[0][:pipeline_finished_at]).to be_nil
    end

    it 'ignores disabled workflows' do
      config, runtime, _ = create_workflow(
        workflow_name: 'my-cat-ingestion',
        run_interval: 0,
        runtime_config:{commit: true},
        run_start: '2025-01-16T12:00:00Z',
        run_stop: '2025-01-16T12:20:00Z',
        disabled: true
      )

      auth_header(:editor)
      get("/api/workflows/labors/status")
      expect(last_response.status).to eq(200)

      expect(json_body.count).to eq(0)
    end

  end
end
