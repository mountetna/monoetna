describe SFTPMetisUploaderJob do
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
        "files_modified_path" => "/tmp/",
        "bucket_name" => "triage",
        "metis_root_path" => "browse/waiting_room/fastq.ucsf.edu"
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
    ]
  }

  let(:fake_sftp_file_stream) {
    StringIO.new("fake file content")
  }

  let(:runtime_config) {
    {
      "commit" => true,
    }
  }

  let(:run_id) { "1234567890" }
  let(:run_record) {
    {
      run_id: run_id,
      config_id: config["config_id"],
      version_number: config["version_number"],
      state: {files_to_update_path: "/tmp/#{run_id}/files_to_update.csv"},
      orchestrator_metadata: {},
      output: nil,
      created_at: Time.now,
      updated_at: Time.now
    }
  }

  def create_files_to_update_csv
    CSV.open(run_record[:state][:files_to_update_path], "wb") do |csv|
      csv << ["path", "modified_time"]
      sftp_files.each do |file|
        csv << [file[:path], file[:modified_time]]
      end
    end
  end
  
  before do
    ENV['TOKEN'] = TEST_TOKEN
    ENV['ARGO_WORKFLOW_ID'] = run_id
  end

  def stub_upload_file_with_stream(file_path, stream)
    stub_upload_file(
      authorize_body: JSON.generate({
        url: "#{METIS_HOST}\/#{config["project_name"]}\/upload/#{config["config"]["metis_root_path"]}/#{file_path}",
      }),
      upload_body: JSON.generate({
        current_byte_position: 0,
        next_blob_size: stream.size,
      }),
      project: config["project_name"],
    )
  end

  context 'records to process' do

    before do
      create_files_to_update_csv
      stub_polyphemus_get_run(config["project_name"], run_id, run_record)

      # sftp stubs
      stub_initial_sftp_connection
      stub_sftp_client_download_as_stream(return_io: fake_sftp_file_stream)
      
      # metis stubs
      stub_metis_setup
      stub_create_folder(bucket: config["config"]["bucket_name"], project: config["project_name"])
      file_to_upload = "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/RA_DB2_SCC1_S32_L007_R1_001.fastq.gz"
      file_to_upload_2 = "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/RA_DB2_SCC1_S32_L007_I1_001.fastq.gz"
      stub_upload_file_with_stream(file_to_upload, fake_sftp_file_stream)
      # TODO: fix the 2nd stub
      stub_upload_file_with_stream(file_to_upload_2, fake_sftp_file_stream)
    end 

    it 'processes the files' do
      job = SFTPMetisUploaderJob.new(config, runtime_config)
      job.execute
    end


  end

end
