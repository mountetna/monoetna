describe Polyphemus::RunJob do
  class DummyJob < Polyphemus::ETLJob
    def pre(context)
      true
    end
    def process(context)
    end
    def post(context)
    end
  end
  it 'runs a job' do
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
           config: {})
    stub_request(:post, "#{JANUS_HOST}/api/tokens/generate").to_return(body: TEST_TOKEN)
    cmd = Polyphemus::RunJob.new
    cmd.setup(Polyphemus.instance.instance_variable_get("@config"))
    cmd.execute('dummy', 1, 1)
  end
end
