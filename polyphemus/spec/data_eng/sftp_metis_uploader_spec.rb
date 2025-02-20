describe SftpMetisUploaderJob do
  include Rack::Test::Methods

  def create_job(config, runtime)
    SftpMetisUploaderJob.new(TEST_TOKEN, config, runtime_config)
  end

  let(:config) {
    config = {
      "project_name" => "labors",
      "secrets" => {
        "sftp_ingest_host" => "some-sftp-host", 
        "sftp_ingest_user" => "user",
        "sftp_ingest_password" => "password"
      },
      "config" => {
        "magic_string" => "LABORS",
        "ingest_root_path" => "SSD",
        "path_to_write_files" => "/tmp/",
        "bucket_name" => "deposit",
        "metis_root_path" => "files/some-sftp-host",
        "deposit_root_path" => "/labors/staging/polyphemus/ingest/some-sftp-host"
      },
      "config_id" => "1",
      "version_number" => "1"
    }
    config
  }

  # TODO: eventually add mutliple files to test with, but running into stubbing and streaming errors
  let(:sftp_files) {
    [
     { path: "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/LABORS_S1.fastq.gz", modified_time: 1672876800 },
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
      state: {
        files_to_update: sftp_files
      },
      orchestrator_metadata: {},
      output: nil,
      created_at: Time.now,
      updated_at: Time.now
    }
  }

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
      file_path: config["config"]["metis_root_path"] + "/" + file_path,
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
    let(:captured_requests) { [] }
    let(:file_to_upload) { "SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/S1.fastq.gz" }
    let(:fake_stream) { @fake_stream ||= StringIO.new("fake file content") }

    before do
      stub_initial_sftp_connection
      stub_polyphemus_get_run(config["project_name"], run_id, run_record)
      stub_sftp_client_download_as_stream(return_io: fake_stream)
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      stub_metis_setup
      stub_create_folder(bucket: config["config"]["bucket_name"], project: config["project_name"])
    end 

    it 'successfully uploads files' do
      stub_upload_file_with_stream(file_to_upload, fake_stream)

      job = create_job(config, runtime_config)
      context = job.execute
      expect(context[:failed_files]).to be_empty
      expect(captured_requests).to be_empty
    end

    it 'fails to upload files' do
      stub_upload_file_with_stream(file_to_upload, fake_stream, force_error: true)

      allow(Polyphemus.instance.logger).to receive(:warn)
      expect(Polyphemus.instance.logger).to receive(:warn).with(/Failed to upload to metis/)
      expect(Polyphemus.instance.logger).to receive(:warn).with(/Found 1 failed files/)
      expect(Polyphemus.instance.logger).to receive(:warn).with(/SSD.*LABORS_S1.fastq.gz/)

      job = create_job(config, runtime_config)
      context = job.execute

      expect(captured_requests[0][:state][:metis_num_failed_files]).to eq(1)
    end
  end

end
