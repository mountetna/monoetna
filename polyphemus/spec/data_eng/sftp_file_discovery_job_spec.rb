describe SftpFileDiscoveryJob do
  include Rack::Test::Methods
  def create_job(config, runtime_config)
    SftpFileDiscoveryJob.new(
      TEST_TOKEN,
      config,
      runtime_config
    )
  end
  let(:config) {
    config = {
      "project_name" => "labors",
      "secrets" => {
        "sftp_ingest_host" => "some-sftp-host", 
        "sftp_ingest_user" => "user",
        "sftp_ingest_password" => "password",
        "sftp_deposit_host" => "other-sftp-host", 
        "sftp_deposit_user" => "user",
        "sftp_deposit_password" => "password",
      },
      "config" => {
        "magic_string" => "LABORS",
        "ingest_root_path" => "SSD"
      },
      "config_id" => "1",
      "version_number" => "1"
    }
  }

  let(:sftp_files) {
    [
      { path: "SSD/ACMK02/LABORS_S1_R1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S1_I1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S1_I2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S1_R2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S2_R1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S2_I1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S2_I2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S2_R2_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S3_R1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S3_I1_001.fastq.gz", modified_time: 1672876800 },
      { path: "SSD/ACMK02/LABORS_S3_I2_001.fastq.gz", modified_time: 1672876800 },
    ]
  }

  let(:run_id) { "1234567890" }

  let(:workflow_name){
    "workflow_name-2000"
  }
  
  before do
      ENV['TOKEN'] = TEST_TOKEN
      ENV['KUBE_ID'] = run_id
      ENV['WORKFLOW_NAME'] = workflow_name
  end

  context 'time window' do

    before do
      stub_initial_sftp_connection
      stub_sftp_search_files(sftp_files)
    end 

    it 'sets the end time to the start time + interval if provided' do
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      runtime_config = {
        "config" => {
          "override_interval" => 864000, # 10 days
          "initial_start_scan_time" => '2023-01-01T00:00:00'
        }
      }
      job = create_job(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:start_time]).to eq(DateTime.parse(runtime_config["config"]["initial_start_scan_time"]).to_i)
      expect(captured_requests[0][:state][:end_time]).to eq(DateTime.parse(runtime_config["config"]["initial_start_scan_time"]).to_i + runtime_config["config"]["override_interval"])
    end

    it 'sets the end time to the current time if no interval is provided' do
      config["config"]["interval"] = nil
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      runtime_config = {
        "config" => {
          "interval" => nil,
          "initial_start_scan_time" => '2023-01-01T00:00:00'
        }
      }
      job = create_job(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:start_time]).to eq(DateTime.parse(runtime_config["config"]["initial_start_scan_time"]).to_i)
      expect(captured_requests[0][:state][:end_time]).to be_within(5).of(Time.now.to_i)
    end

  end


  context 'first run (restarted)' do

    let(:runtime_config) {
      {
        "config" => {
          "interval" => 864000, # 10 days
          "initial_start_scan_time" => '2023-01-01T00:00:00'
        }
      }
    }

    before do
      stub_initial_sftp_connection
      stub_sftp_search_files(sftp_files)
    end

    it "sends the proper run state to the db" do
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      job = create_job(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:files_to_update].size).to eq(11)
      expect(captured_requests[0][:state][:files_to_update]).to all(include(:path))
    end
  end

  context 'subsequent runs' do


    let(:last_state) { {"end_time" => 1672876800} } # Jan 5, 2023
    let(:runtime_config) {
      {
        "config" => {
          "override_interval" => 864000, # 10 days
          "initial_start_scan_time" => nil, #Jan 1, 2023
        }
      }
    }

    before do
      stub_polyphemus_get_last_state(
        config["project_name"],
        config["config_id"],
        last_state
      )
      stub_initial_sftp_connection
      stub_sftp_search_files(sftp_files)
    end

    it "fetches the last scan timestamp from the database" do
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      job = create_job(config, runtime_config)
      job.execute
      expect(captured_requests[0][:state][:start_time]).to eq(last_state["end_time"])
      expect(captured_requests[0][:state][:end_time]).to eq(last_state["end_time"] + 864000)
    end
  end
end
