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
        "path_to_write_files" => "/tmp/",
        "bucket_name" => "triage",
        "metis_root_path" => "browse/waiting_room/fastq.ucsf.edu"
      },
      "config_id" => "1",
      "version_number" => "1"
    }
    config
  }

  # TODO: eventually add mutliple files to test with, but running into stubbing and streaming errors
  let(:sftp_files) {
    [
     { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/DSCOLAB_RA_DB2_SCC1_S32_L007_R1_001.fastq.gz", modified_time: 1672876800 },
    ]
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
      state: {files_to_update_path: "/tmp/#{run_id}/#{SFTPFileDiscoveryJob::SFTP_FILES_TO_UPDATE_CSV}"},
      orchestrator_metadata: {},
      output: nil,
      created_at: Time.now,
      updated_at: Time.now
    }
  }

  def create_files_to_update_csv
    Dir.mktmpdir(run_id, "/tmp/#{run_id}")
    CSV.open(run_record[:state][:files_to_update_path], "wb") do |csv|
      csv << ["path", "modified_time"]
      sftp_files.each do |file|
        csv << [file[:path], file[:modified_time]]
      end
    end
  end


  let(:workflow_name){
    "workflow_name-2000"
  }
  
  before do
      ENV['TOKEN'] = TEST_TOKEN
      ENV['KUBE_ID'] = run_id
      ENV['WORKFLOW_NAME'] = workflow_name
  end

  def stub_upload_file_with_stream(file_path, stream, force_error: false)
    next_blob_size = force_error ? stream.size * 100000: stream.size
    stub_upload_file(
      authorize_body: JSON.generate({
        url: "#{METIS_HOST}\/#{config["project_name"]}\/upload/#{config["config"]["metis_root_path"]}/#{file_path}",
      }),
      upload_body: JSON.generate({
        current_byte_position: 0,
        next_blob_size: next_blob_size,
      }),
      project: config["project_name"],
    )
  end

  context 'records to process' do

    before do
      create_files_to_update_csv
      stub_initial_sftp_connection
      stub_polyphemus_get_run(config["project_name"], run_id, run_record)
    end 

    it 'successfully uploads files' do
      fake_stream = StringIO.new("fake file content")
      stub_sftp_client_download_as_stream(return_io: fake_stream)

      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)

      # metis stubs
      stub_metis_setup
      stub_create_folder(bucket: config["config"]["bucket_name"], project: config["project_name"])
      # these files have the DSCOLAB prefix removed
      file_to_upload = "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/RA_DB2_SCC1_S32_L007_R1_001.fastq.gz"
      stub_upload_file_with_stream(file_to_upload, fake_stream)

      job = SFTPMetisUploaderJob.new(config, runtime_config)
      context = job.execute
      expect(context[:failed_files]).to be_empty
      expect(captured_requests).to be_empty
    end

    it 'fails to upload files' do
      fake_stream = StringIO.new("fake file content")
      stub_sftp_client_download_as_stream(return_io: fake_stream)
      captured_requests = []
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)

      # metis stubs
      stub_metis_setup
      stub_create_folder(bucket: config["config"]["bucket_name"], project: config["project_name"])
      file_to_upload = "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/RA_DB2_SCC1_S32_L007_R1_001.fastq.gz"
      # here we just force an streaming error to trigger an exception
      stub_upload_file_with_stream(file_to_upload, fake_stream, force_error: true)

      job = SFTPMetisUploaderJob.new(config, runtime_config)
      context = job.execute
      failed_files_path = "/tmp/1234567890/#{SFTPMetisUploaderJob::METIS_FAILED_FILES_CSV}"
      expect(File.exist?(failed_files_path)).to be_truthy
      expect(captured_requests[0][:state][:metis_num_failed_files]).to eq(1)
      expect(captured_requests[0][:state][:metis_failed_files_path]).to eq(failed_files_path)
    end
  end

end