describe Polyphemus::RuntimeConfig do

  let(:some_config) do
    {
      "some_key" => "some_value"
    }
  end

  before do
    Timecop.freeze('2025-01-16T12:30:00Z')
    create_workflow(
      workflow_name: 'ready-workflow',
      run_interval: 1000,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: '2025-01-16T12:10:00Z',
    )

    create_workflow(
      workflow_name: 'running-workflow',
      run_interval: 1000,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: nil,
      status: 'Running'
    )

    create_workflow(
      workflow_name: 'pending-workflow',
      run_interval: 1000_000,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: '2025-01-16T12:10:00Z',
    )

    create_workflow(
      workflow_name: 'failed-workflow',
      run_interval: 1000,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: '2025-01-16T12:20:00Z',
      status: "Failed"
    )

    fixed_config, *_ = create_workflow(
      workflow_name: 'fixed-workflow',
      run_interval: 1000,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: '2025-01-16T12:05:00Z',
      status: "Failed"
    )
    create_run(fixed_config, 1001,
      run_start: '2025-01-16T12:06:00Z',
      run_stop: '2025-01-16T12:08:00Z',
      status: 'Succeeded'
    )

    create_workflow(
      workflow_name: 'ignored-workflow',
      run_interval: 0,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: '2025-01-16T12:10:00Z',
    )

    create_workflow(
      workflow_name: 'ignored-workflow-2',
      run_interval: nil,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: '2025-01-16T12:10:00Z',
    )

    create_workflow(
      workflow_name: 'disabled-workflow',
      run_interval: 100000,
      disabled: true,
      run_start: '2025-01-16T12:00:00Z',
      run_stop: '2025-01-16T12:10:00Z',
    )

    create_workflow(
      workflow_name: 'new-workflow',
      run_interval: 1000
    )
  end

  after do
    Timecop.return
  end

  it 'correctly returns configs eligible for scheduling' do
    eligible_runtime_configs = Polyphemus::RuntimeConfig.eligible_runtime_configs
    expect(eligible_runtime_configs.length).to eq(3)
    expect(eligible_runtime_configs.map(&:workflow_config).map(&:workflow_name)).to match_array(
      ["ready-workflow", "new-workflow", "fixed-workflow"]
    )
  end
end
