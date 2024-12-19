describe SFTPFileDiscoveryJob do
  include Rack::Test::Methods
  let(:config) {
    config = {
      "project_name" => "labors",
      "secrets" => {
        "sftp_host" => "monoetna-polyphemus_sftp-1", 
        "sftp_user" => "user",
        "sftp_password" => "password",
        "sftp_port" => "22"
      },
      "config" => {
        "regex" => ".*",
        "root_dir" => "/home/user/uploads"
      },
      "config_id" => "1",
      "version_number" => "1"
    }
    config
  }

  let(:runtime_config) {
    {
      "commit" => true,
    }
  }

  let(:run_id) { "1234567890" }
  let(:empty_run_record) {
    Polyphemus::Run.new(
      run_id: run_id,
      config_id: config["config_id"],
      version_number: config["version_number"],
      state: {},
      orchestrator_metadata: {},
      output: nil,
      created_at: Time.now,
      updated_at: Time.now
    )
  }



  before do
    ENV['TOKEN'] = TEST_TOKEN
    ENV['ARGO_WORKFLOW_ID'] = run_id
    # Setup sftp client
  end

  context 'first run' do

    before do
      WebMock.enable!
      require 'pry'; binding.pry
      stub_polyphemus_get_previous_run(config["project_name"], config["config_id"], config["version_number"], empty_run_record)
      job = SFTPFileDiscoveryJob.new(config, runtime_config)
      job.process 
    end

    it "records the last scan timestamp in the db" do
      # check the db
    end

    it 'records the number of files to update in the db' do
    end

    it 'writes the files to update to a csv' do
    end 


  end

  context 'subsequent runs' do

    before do
      # set the db
    end

    it "fetches the last scan timestamp from the database" do
      job = SFTPFileDiscoveryJob.new(config, secrets)
      job.fetch_last_scan
    end
  end

end
