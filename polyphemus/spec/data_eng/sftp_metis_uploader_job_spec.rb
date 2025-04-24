describe SftpMetisUploaderJob do
  include Rack::Test::Methods

  def create_job(config, runtime)
    SftpMetisUploaderJob.new(TEST_TOKEN, config, runtime)
  end

  let(:config) {
    config = {
      "project_name" => "labors",
      "secrets" => {
        "sftp_ingest_host" => "some-sftp-host", 
        "sftp_ingest_user" => "user",
        "sftp_ingest_password" => "password",
        "notification_webhook_url" => "https://deposit.slack.test/DEPOSIT/DEPOSIT"
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
      'config' => {
        "commit" => true,
      }
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
    let(:fake_stream) { @fake_stream ||= StringIO.new("A"*32) }

    before do
      stub_initial_sftp_connection
      stub_polyphemus_get_run(config["project_name"], run_id, run_record)
      stub_sftp_client_download_as_stream(return_io: fake_stream)
      stub_polyphemus_update_run(config["project_name"], run_id, captured_requests)
      stub_metis_setup
      stub_create_folder(bucket: config["config"]["bucket_name"], project: config["project_name"])
      
      stub_request(
        :get, "http://sftp//some-sftp-host/SSD/20240919_LH00416_0184_B22NF2WLT3/ACMK02/"
      ).to_return(
        body: sftp_files.map { |f|
          "drwxrwxr-x    4 eurystheus     labors        32 Aug 26  2022 #{::File.basename(f[:path])}"
        }.join("\n"),
        headers: {}
      )
      stub_slack(config["secrets"]["notification_webhook_url"])
      stub_polyphemus_get_last_state(
        config["project_name"],
        config["config_id"],
        {}
      )
    end

    it 'successfully uploads files' do
      stub_upload_file_with_stream(file_to_upload, fake_stream)
      allow_any_instance_of(Etna::Filesystem::SftpFilesystem).to receive(:with_readable).and_yield(StringIO.new("A"*32))

      job = create_job(config, runtime_config)
      context = job.execute
      expect(context[:failed_files]).to be_empty
      expect(captured_requests.first[:state]).to eq(
        metis_successful_files: [ sftp_files.first[:path] ]
      )
    end

    it 'skips already uploaded files' do
      stub_upload_file_with_stream(file_to_upload, fake_stream)
      allow_any_instance_of(Etna::Filesystem::SftpFilesystem).to receive(:with_readable).and_yield(StringIO.new("A"*32))

      stub_polyphemus_get_last_state(
        config["project_name"],
        config["config_id"],
        {
          metis_successful_files: [ [ sftp_files.first[:path] ] ]
        }
      )

      job = create_job(config, runtime_config)
      context = job.execute

      expect(context[:files_to_update]).to be_empty
      expect(captured_requests).to eq([])
    end

    context 'override root path' do
      let(:override_root_path) { "archive" }
      let(:alt_runtime_config) {
        {
          'config' => { 'commit' => true, 'override_root_path' => override_root_path } 
        }
      }
      let(:alt_run_record) {
        run_record.merge(
          state: {
            files_to_update: sftp_files.map do |file|
              file.merge(
                path: file[:path].sub(/^SSD/,override_root_path)
              )
            end
          }
        )
      }

      it 'skips already uploaded files in the original root path' do
        stub_polyphemus_get_run(config["project_name"], run_id, alt_run_record)
        allow_any_instance_of(Etna::Filesystem::SftpFilesystem).to receive(:with_readable).and_yield(StringIO.new("A"*32))

        stub_polyphemus_get_last_state(
          config["project_name"],
          config["config_id"],
          {
            metis_successful_files: [ [ sftp_files.first[:path] ] ]
          }
        )

        job = create_job(config, alt_runtime_config)
        context = job.execute

        expect(context[:files_to_update]).to be_empty
        expect(captured_requests).to eq([])
      end
    end

    it 'fails to upload files' do
      allow_any_instance_of(Etna::Filesystem::SftpFilesystem).to receive(:with_readable).and_yield(StringIO.new)
      stub_upload_file_with_stream(file_to_upload, fake_stream)

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
