describe WorkflowController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  before(:each) do
    stub_event_log
  end

  let(:run_id) { "argo_run_id_22" }

  def create_workflow_with_configs
    auth_header(:editor)
    json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
    expect(last_response.status).to eq(200)
    some_config = { "test_key" => 'update 1' }
    json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
      workflow_name: 'my workflow name',
      workflow_type: 'test-workflow',
      config: some_config
    )
    expect(last_response.status).to eq(200)
  end

  context 'config creation' do
    it 'initializes a workflow' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')

      expect(last_response.status).to eq(200)
      expect(Polyphemus::Config.count).to eq(1)
      config = Polyphemus::Config.last

      expect(config.workflow_name).to eq('my workflow name')
      expect(config.workflow_type).to eq('test-workflow')
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
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
      expect(last_response.status).to eq(200)
      expect(Polyphemus::RuntimeConfig.count).to eq(1)
      expect(Polyphemus::RuntimeConfig.last.config_id).to eq(json_body[:config_id])
    end
  end

  context 'config updates' do

    it 'correctly updates the config and version number' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'test_key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        config: some_config
      )
      expect(last_response.status).to eq(200)
      first_config = Polyphemus::Config.last
      # Contains the updated config
      expect(first_config.config).to eq(some_config)
      # Contains the same config_id
      expect(first_config.config_id).to eq(json_body[:config_id])
      # Contains a new version number
      expect(first_config.version_number).to eq(2)


      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        config: another_config
      )
      expect(last_response.status).to eq(200)
      second_config = Polyphemus::Config.last
      expect(second_config.config_id).to eq(json_body[:config_id])
      expect(second_config.config).to eq(another_config)
      expect(second_config.version_number).to eq(3)
      expect(second_config.updated_at).to_not eq(first_config.updated_at)
      expect(second_config.created_at).to_not eq(first_config.created_at)
    end

    it 'rejects an invalid config' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'not a valid key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        config: some_config
      )
      expect(last_response.status).to eq(422)
    end

    it 'updates secrets and not the version number' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
      expect(last_response.status).to eq(200)
      some_secret = { 'test_secret' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        secrets: some_secret
      )
      expect(last_response.status).to eq(200)
      the_config = Polyphemus::Config.last
      expect(the_config.secrets).to eq(some_secret)
      expect(the_config.version_number).to eq(1)

      another_secret = { 'test_secret' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        secrets: another_secret
      )
      expect(last_response.status).to eq(200)
      updated_config = Polyphemus::Config.last
      expect(updated_config.config_id).to eq(json_body[:config_id])
      expect(updated_config.secrets).to eq(another_secret)
      expect(updated_config.version_number).to eq(1)
      expect(updated_config.created_at).to eq(the_config.created_at)
      expect(updated_config.updated_at).to_not eq(the_config.updated_at)
    end

    it 'rejects unknown secrets' do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
      expect(last_response.status).to eq(200)
      some_secret = { 'not_a_secret' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        secrets: some_secret
      )
      config = Polyphemus::Config.last
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq('Secrets for test-workflow must be one of: test_secret')
      expect(config.secrets).to eq({})
    end

  end

  context 'retrieve a config' do
    before do
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'test_key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        config: some_config
      )
      expect(last_response.status).to eq(200)
      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        config: another_config
      )
      expect(last_response.status).to eq(200)
    end

    it 'returns the most recent config when version is not specified' do
      auth_header(:editor)
      get("/api/workflows/labors/configs/#{json_body[:config_id]}")
      expect(last_response.status).to eq(200)
      expect(json_body[:version_number]).to eq(3)
    end

    it 'returns a config by id and version' do
      auth_header(:editor)
      get("/api/workflows/labors/configs/#{json_body[:config_id]}?version=1")
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
        workflow_type: 'test-workflow',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )
      config = Polyphemus::Config.create(
        project_name: 'athena',
        workflow_name: 'my 2nd workflow name',
        workflow_type: 'test-workflow',
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
        workflow_type: 'test-workflow',
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
      post('/api/workflows/configs', workflow_type: 'test-workflow')
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
      auth_header(:editor)
      json_post('/api/workflows/labors/create', workflow_name: 'my workflow name', workflow_type: 'test-workflow')
      expect(last_response.status).to eq(200)
      some_config = { 'test_key' => 'update 1' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
        config: some_config
      )
      expect(last_response.status).to eq(200)
      another_config = { 'test_key' => 'update 2' }
      json_post("/api/workflows/labors/update/#{json_body[:config_id]}",
        workflow_name: 'my workflow name',
        workflow_type: 'test-workflow',
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
        config_id: json_body[:config_id],
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
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        name: "test-workflow-1999",
        output: initial_output
      )
      expect(last_response.status).to eq(200)

      # Update with additional output
      json_post("/api/workflows/labors/run/update/#{run_id}",
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
    before do
      create_workflow_with_configs
      # Create a new run
      auth_header(:editor)
      the_state = {:some_state => "ya_hooo!", :another_state => "woohoo!"}
      json_post("/api/workflows/labors/run/update/#{run_id}",
        config_id: json_body[:config_id],
        version_number: 2,
        name: "test-workflow-1999",
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

    it 'updates existing runtime config' do
      config = {:commit  => true }
      json_post("/api/workflows/labors/runtime_configs/update/#{json_body[:config_id]}",
        config: config,
        run_interval: 60, 
        disabled: true
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
      expect(json_body[:disabled]).to eq(true)
    end

    it 'raises an error when the runtime config is invalid' do
      config = {:not_a_valid_param  => "not a boolean" }
      json_post("/api/workflows/labors/runtime_configs/update/#{json_body[:config_id]}",
        config: config,
      )
      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq("no such param not_a_valid_param")
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
      expect(json_body[:error]).to eq('No such config 10000') 
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

  context 'running a workflow once' do
    before do
      create_workflow_with_configs
    end

    it 'successfully runs a workflow' do
      config = Polyphemus::Config.current.where(project_name: 'labors', config_id: json_body[:config_id]).first
      
      # Stub argo client command
      workflow_output = "Workflow Name: simple-two-jobs-ktm4g\nWorkflow UID: de6df0ed-b32c-4707-814f-11c323b0687b\n"
      allow(Open3).to receive(:capture3).and_return([workflow_output, "", double(success?: true)])
      json_post("/api/workflows/labors/runtime_configs/run_once/#{config.config_id}", 
        workflow_type: "test-workflow",
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
      config = Polyphemus::Config.current.where(project_name: 'labors', config_id: json_body[:config_id]).first
      allow(Open3).to receive(:capture3).and_return(["", "Command failed", double(success?: false)])
      json_post("/api/workflows/labors/runtime_configs/run_once/#{config.config_id}", 
        workflow_type: "test-workflow",
        config: { commit: true }
      )

      expect(last_response.status).to eq(422)
      expect(json_body[:error]).to eq("Failed to submit Argo workflow: Command failed")
    end

  end

  context 'status' do

    it 'returns the status of a workflow that has completed' do

      config_id = Polyphemus::Config.next_id

      config_1 = Polyphemus::Config.create(
        project_name: 'labors',
        workflow_name: 'my-cat-ingestion',
        workflow_type: 'cat-ingestion',
        config_id: config_id,
        version_number: 1,
        config: {"some_key" => "some_value"},
        secrets: {},
        created_at: '2025-01-16T12:00:00Z',
        updated_at: '2025-01-16T12:00:00Z'
      )

      config_2 = Polyphemus::Config.create(
        project_name: 'labors',
        workflow_name: 'my-cat-ingestion',
        workflow_type: 'cat-ingestion',
        config_id: config_id,
        version_number: 2,
        config: {"some_key" => "some_value_2"},
        secrets: {},
        created_at: '2025-01-18T12:00:00Z',
        updated_at: '2025-01-18T12:00:00Z'
      )

      runtime_config_1 = Polyphemus::RuntimeConfig.create(
        config_id: config_id,
        run_interval: 1000,
        config:{commit: true} 
      )

      # Create a run for the first version
      Polyphemus::Run.create(
        run_id: "hello-there-1",
        config_id: config_id,
        version_number: config_1.version_number,
        name: 'my-cat-ingestion-1000',
        orchestrator_metadata: {
          'startedAt' => '2025-01-16T12:00:00Z',
          'nodes' => {
            'steps-node' => {
              'finishedAt' => '2025-01-16T12:20:00Z',
              'phase' => 'Succeeded',
              'type' => 'Steps'
            }
          }
        }
      )

      # Create a run for the second version
      Polyphemus::Run.create(
        run_id: "hello-there-2",
        config_id: config_id,
        version_number: config_2.version_number,
        name: 'my-cat-ingestion-2000',
        orchestrator_metadata: {
          'startedAt' => '2025-01-18T12:00:00Z',
          'nodes' => {
            'steps-node' => {
              'finishedAt' => '2025-01-18T12:20:00Z',
              'phase' => 'Succeeded',
              'type' => 'Steps'
            }
          }
        }
      )

      auth_header(:editor)
      get("/api/workflows/labors/status")
      expect(last_response.status).to eq(200)

      # Make sure the second version of the config is used (we use the latest run)
      expect(json_body.count).to eq(1)
      expect(json_body[0][:workflow_name]).to eq("my-cat-ingestion")
      expect(json_body[0][:workflow_type]).to eq("cat-ingestion")
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

    it 'returns a nil status when the workflow hasnt been run yet' do
      config_2 = Polyphemus::Config.create(
        project_name: 'labors',
        workflow_name: 'my-cat-ingestion-2',
        workflow_type: 'cat-ingestion',
        config_id: Polyphemus::Config.next_id,
        version_number: 1,
        config: {},
        secrets: {},
      )

      auth_header(:editor)
      get("/api/workflows/labors/status")
      expect(last_response.status).to eq(200)
      expect(json_body.count).to eq(1)
      expect(json_body[0][:workflow_name]).to eq("my-cat-ingestion-2")
      expect(json_body[0][:workflow_type]).to eq("cat-ingestion")
      expect(json_body[0][:pipeline_state]).to be_nil
      expect(json_body[0][:pipeline_finished_at]).to be_nil
    end
  end
end
