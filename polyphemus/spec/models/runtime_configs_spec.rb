describe Polyphemus::RuntimeConfig do

  let(:some_config) do
    {
      "some_key" => "some_value"
    }
  end

  before do
    config_1 = Polyphemus::Config.create(
      project_name: 'labors',
      workflow_name: 'my-cat-ingestion',
      workflow_type: 'cat-ingestion',
      config_id: Polyphemus::Config.next_id,
      version_number: 1,
      config: {},
      secrets: {},
    )
    runtime_config_1 = Polyphemus::RuntimeConfig.create(
      config_id: config_1.config_id,
      run_interval: 1000,
      config: some_config
    )

    Polyphemus::Run.create(
      run_id: "hello-there-1",
      config_id: config_1.config_id,
      version_number: config_1.version_number,
      name: 'my-cat-ingestion-1000',
      orchestrator_metadata: {
        'startedAt' => '2025-01-16T12:00:00Z',
        "nodes" => {
          "step-node" => {
            'finishedAt' => '2025-01-16T12:20:00Z',
            'type' => 'Steps',
            'phase' => 'Succeeded'
          }
        }
      }
    )

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

    Polyphemus::Run.create(
      run_id: "hello-there-2",
      config_id: config_2.config_id,
      version_number: config_2.version_number,
      name: 'my-cat-ingestion-2000',
      orchestrator_metadata: {
        'startedAt' => '2025-01-16T12:00:00Z',
        "nodes" => {
          "step-node" => {
            'finishedAt' => '2025-01-16T12:10:00Z',
            'type' => 'Steps',
            'phase' => 'Succeeded'
          }
        }
      }
    )

  end

  it 'correctly returns configs eligible for scheduling' do
      eligible_runtime_configs = Polyphemus::RuntimeConfig.eligible_runtime_configs
      expect(eligible_runtime_configs[0].config).to eq(some_config)
  end
end
