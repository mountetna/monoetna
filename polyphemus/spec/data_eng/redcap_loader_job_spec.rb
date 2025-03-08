describe RedcapLoaderJob do
  def create_job(config, runtime_config)
    RedcapLoaderJob.new(
      TEST_TOKEN,
      config,
      runtime_config
    )
  end

  context 'loads' do
    let(:config) do
      {
        'project_name' => 'test',
        'config' => {},
        'secrets' => {
          'redcap_tokens' => REDCAP_TOKEN
        }
      }
    end
    let(:runtime_config) do
      {
        'config' => {
          'mode' => 'append',
          'model_names' => 'all'
        }
      }
    end
    let(:run_id) { "1234567890" }

    it 'all models' do
      ENV['KUBE_ID'] = run_id
      ENV['TOKEN'] = TEST_TOKEN
      stub_magma_models
      copy_redcap_project
      stub_redcap_data(:essential_data)
      stub_request(:post, "#{POLYPHEMUS_HOST}/api/workflows/test/run/update/#{run_id}").to_return(body: "{}")
      job = create_job( config, runtime_config )

      expect {
        context = job.execute
      }.not_to raise_error
    end
  end
end
