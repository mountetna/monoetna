describe Polyphemus::IntervalScheduler do
  class DummyManifest < Polyphemus::WorkflowManifest
    class << self
      def as_json
        {
          name: 'dummy',
          schema: {
            "$schema": "http://json-schema.org/draft-07/schema#",
            type: "object",
            description: "Just a dummy",
	    title: "Dummy Loader"
          },
          runtime_params: { },
          secrets: [ ],
          workflow_path: '/app/workflows/argo/dummy/workflow.yaml'
        }
      end
    end
  end
  class DummyJob < Polyphemus::ETLJob
    def pre(context)
      true
    end
    def process(context)
    end
    def post(context)
    end
  end
  it 'finds scheduled jobs' do
    create(:config,
           project_name: 'labors',
           workflow_name: 'dummy workflow',
           workflow_type: 'dummy',
           version_number: 1,
           config_id: 1,
           secrets: {},
           config: {})
    create(:runtime_config,
           config_id: 1,
           run_interval: 1000,
           config: {})
    cmd = Polyphemus::IntervalScheduler.new
    cmd.setup(Polyphemus.instance.instance_variable_get("@config"))

    stub_const("Polyphemus::IntervalScheduler::SLEEP_INTERVAL", 0.1)
    allow_any_instance_of(Polyphemus::IntervalScheduler).to receive(:loop).and_yield
    allow(Polyphemus.instance.logger).to receive(:info)
    expect(Polyphemus.instance.logger).not_to receive(:error)

    status = double()
    allow(status).to receive(:success?).and_return(true)
    allow(Open3).to receive(:capture3).and_return(
      [
        "Workflow Name: dummy\nWorkflow UID: dummy_id",
        "",
        status
      ]
    )

    cmd.execute

    expect(Polyphemus.instance.logger).to have_received(:info).with(/Found 1 eligible runtime configs for scheduling/)
    expect(Polyphemus.instance.logger).to have_received(:info).with(/Submitting workflow dummy workflow, for project: labors, workflow_type: dummy, config_id: 1/)
    expect(Polyphemus.instance.logger).to have_received(:info).with(/Submitting Argo workflow with command:/)
    expect(Polyphemus.instance.logger).to have_received(:info).with(/Argo workflow submitted successfully - Name: dummy, UID: dummy_id/)
    expect(Polyphemus::Run.count).to eq(1)
  end
end
