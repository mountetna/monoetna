describe Etna::Clients::Metis::IngestMetisDataWorkflow do
  describe "with SFTP" do
    let(:ingest_filesystem) {
      Etna::Filesystem::Mock.new do |fname, opts|
        StringIO.new
      end
    }
    let(:metis_client) {
      Etna::Clients::Metis.new(
        token: TEST_TOKEN,
        host: METIS_HOST,
      )
    }
    let(:metis_filesystem) {
      Etna::Filesystem::Metis.new(
        project_name: PROJECT,
        bucket_name: "triage",
        metis_client: metis_client,
      )
    }

    before(:each) do
      file_contents = "simple-file-contents"
      stub_metis_setup
      stub_create_folder(bucket_name: "waiting_room")
      stub_upload_file(
        authorize_body: JSON.generate({
          url: "#{METIS_HOST}\/#{PROJECT}\/upload/test.txt",
        }),
        upload_body: JSON.generate({
          current_byte_position: file_contents.bytesize,
          next_blob_size: 0,
        }),
      )

      ingest_filesystem.with_writeable("test.txt") do |f|
        f.write(file_contents)
      end
    end

    it "can ingest to Metis" do
      workflow = Etna::Clients::Metis::IngestMetisDataWorkflow.new(
        metis_filesystem: metis_filesystem,
        ingest_filesystem: ingest_filesystem,
        logger: nil,
      )

      workflow.copy_files(["test.txt"])

      expect(WebMock).to have_requested(:post, "#{METIS_HOST}/authorize/upload")
      # Once to start the upload, once to send the blob.
      expect(WebMock).to have_requested(:post, "#{METIS_HOST}/#{PROJECT}/upload/test.txt").times(2)
    end
  end
end
