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
        "path_to_write_files" => "/tmp/" 
      },
      "config_id" => "1",
      "version_number" => "1"
    }
    config
  }

  let(:last_scan_timestamp) {
    "1672531200" #Jan 1, 2023
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
  let(:empty_run_record) {
    {}
  }
  let(:run_record) {
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
  end

  context 'first run' do

    before do
      stub_polyphemus_get_previous_run(config["project_name"], config["config_id"], config["version_number"], empty_run_record)
      stub_initial_sftp_connection
      stub_sftp_search_files(sftp_files)
      job = SFTPFileDiscoveryJob.new(config, runtime_config)
      job.execute
    end

    it "records the last scan timestamp in the db" do
    end

    it 'records the number of files to update in the db' do
    end

    it 'writes the files to update to a csv' do
    end 


  end

  context 'subsequent runs' do

    before do
    end

    it "fetches the last scan timestamp from the database" do
    end
  end

end
