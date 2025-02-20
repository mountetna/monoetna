describe SftpDepositUploaderJob do
    include Rack::Test::Methods
    def create_job(config, runtime_config)
      SftpDepositUploaderJob.new(TEST_TOKEN, config, runtime_config)
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
          "sftp_deposit_password" => "password"
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

    let(:workflow_name){
      "workflow_name-2000"
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
  
    before do
      ENV['TOKEN'] = TEST_TOKEN
      ENV['KUBE_ID'] = run_id
      ENV['WORKFLOW_NAME'] = workflow_name
    end
  
    context 'records to process' do
  
      let(:captured_requests) { [] }

      before do
        stub_initial_sftp_connection
        stub_initial_ssh_connection
        stub_polyphemus_get_run(config["project_name"], run_id, run_record)
        fake_stream = StringIO.new("fake file content")
        stub_sftp_client_download_as_stream(return_io: fake_stream)
        stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      end 
  
      it 'successfully uploads files' do
        stub_remote_ssh_file_upload(success: true)
  
        stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
  
        job = create_job(config, runtime_config)
        context = job.execute

        expect(context[:failed_files]).to be_empty
        expect(captured_requests).to be_empty
      end
  
      it 'fails to upload files' do
        allow(Polyphemus.instance.logger).to receive(:warn)
        expect(Polyphemus.instance.logger).to receive(:warn).with(/Failed to upload to deposit host other-sftp-host/)
        expect(Polyphemus.instance.logger).to receive(:warn).with(/Found 1 failed files/)
        expect(Polyphemus.instance.logger).to receive(:warn).with(/SSD.*LABORS_S1.fastq.gz/)
        stub_remote_ssh_file_upload(success: false)

        job = create_job(config, runtime_config)
        context = job.execute

        expect(captured_requests[0][:state][:deposit_num_failed_files]).to eq(1)
      end
  
    end
  
  end
