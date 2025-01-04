describe SFTPFileDiscoveryJob do
  include Rack::Test::Methods
  let(:config) {
    config = {
      "project_name" => "labors",
      "secrets" => {
        "sftp_host" => "some-sftp-host", 
        "sftp_user" => "user",
        "sftp_password" => "password",
        "sftp_port" => "22"
      },
      "config" => {
        "file_regex" => "DSCOLAB(-|_).*",
        "sftp_root_dir" => "SSD",
        "interval" => 864000, # 10 days
        "initial_start_scan_time" => 1672531200, #Jan 1, 2023
        "files_modified_path" => "/tmp/" 
      },
      "config_id" => "1",
      "version_number" => "1"
    }
    config
  }

  let(:sftp_files) {
    [
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB2_SCC1_S32_L007_R1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB2_SCC1_S32_L007_I1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB2_SCC1_S32_L007_I2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB2_SCC1_S32_L007_R2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_IPIMEL575-T1-CD45enriched-SCC1_S33_L007_R1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_IPIMEL575-T1-CD45enriched-SCC1_S33_L007_I1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_IPIMEL575-T1-CD45enriched-SCC1_S33_L007_I2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_IPIMEL575-T1-CD45enriched-SCC1_S33_L007_R2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB1_SCG1_S22_L008_R1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB1_SCG1_S22_L008_I1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB1_SCG1_S22_L008_I2_001.fastq.gz", modified_time: 1672876800 },
    ]
  }

  let(:runtime_config) {
    {
      "commit" => true,
    }
  }
  let(:run_id) { "1234567890" }
  let(:empty_run_record) {
    {}
  }
  let(:run_record) {
    {
      run_id: run_id,
      config_id: config["config_id"],
      version_number: config["version_number"],
      state: {end_time: 1672876800}, # Jan 5, 2023
      orchestrator_metadata: {},
      output: nil,
      created_at: Time.now,
      updated_at: Time.now
    }
  }
  
  before do
    ENV['TOKEN'] = TEST_TOKEN
    ENV['ARGO_WORKFLOW_ID'] = run_id
  end

  context 'time window' do

    before do
      stub_polyphemus_get_previous_run(config["project_name"], config["config_id"], config["version_number"], empty_run_record)
      stub_initial_sftp_connection
      stub_sftp_search_files(sftp_files)
    end 

    it 'sets the end time to the start time + interval if provided' do
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      job = SFTPFileDiscoveryJob.new(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:start_time]).to eq(config["config"]["initial_start_scan_time"])
      expect(captured_requests[0][:state][:end_time]).to eq(config["config"]["initial_start_scan_time"] + config["config"]["interval"])
    end

    it 'sets the end time to the current time if no interval is provided' do
      config["config"]["interval"] = nil
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      job = SFTPFileDiscoveryJob.new(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:start_time]).to eq(config["config"]["initial_start_scan_time"])
      expect(captured_requests[0][:state][:end_time]).to be_within(5).of(Time.now.to_i)
    end

  end


  context 'first run' do

    before do
      stub_polyphemus_get_previous_run(config["project_name"], config["config_id"], config["version_number"], empty_run_record)
      stub_initial_sftp_connection
      stub_sftp_search_files(sftp_files)
    end

    it "sends the proper run state to the db" do
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      job = SFTPFileDiscoveryJob.new(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:num_files_to_update]).to eq(11)
      expect(captured_requests[0][:state][:files_to_update_path]).to eq("/tmp/1234567890-files_to_update.txt")
    end

    it 'writes the files to update to a csv' do
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      job = SFTPFileDiscoveryJob.new(config, runtime_config)
      job.execute
      expect(File.exist?(captured_requests[0][:state][:files_to_update_path])).to be_truthy
      expect(File.read(captured_requests[0][:state][:files_to_update_path])).to eq("path,modified_time\n" + sftp_files.map { |file| "#{file[:path]},#{file[:modified_time]}" }.join("\n") + "\n")
    end 

  end

  context 'subsequent runs' do

    before do
      stub_polyphemus_get_previous_run(config["project_name"], config["config_id"], config["version_number"], run_record)
      stub_initial_sftp_connection
      stub_sftp_search_files(sftp_files)
    end

    it "fetches the last scan timestamp from the database" do
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      job = SFTPFileDiscoveryJob.new(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:start_time]).to eq(run_record[:state][:end_time])
      expect(captured_requests[0][:state][:end_time]).to eq(run_record[:state][:end_time] + 864000)
    end
  end

end
