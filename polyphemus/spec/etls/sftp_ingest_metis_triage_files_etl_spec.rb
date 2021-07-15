describe Polyphemus::SftpIngestMetisTriageFilesEtl do
  let(:file_contents) {
    "testing, 1, 2, 3!"
  }

  before(:each) do
    stub_data
    stub_metis_setup
    stub_create_folder({
      project: "triage",
      bucket: "waiting_room",
    })
    stub_upload_file({
      project: "triage",
      authorize_body: JSON.generate({
        url: "#{METIS_HOST}\/triage\/upload/foo/bar/test3.txt",
      }),
      upload_body: JSON.generate({
        current_byte_position: file_contents.bytesize,
        next_blob_size: 0,
      }),
    })
  end

  it "should send files to Metis and update timestamp" do
    stub_ingest_filesystem

    ingest_etl = Polyphemus::SftpIngestMetisTriageFilesEtl.new

    expect(Polyphemus::IngestFile.find(name: /test3.txt/)[:ingested_at]).to eq(nil)

    ingest_etl.run_once

    # Only one file should be uploaded
    expect(WebMock).to have_requested(:post, "#{METIS_HOST}/authorize/upload")
    # Once to start the upload, once to send the blob.
    expect(WebMock).to have_requested(:post, "#{METIS_HOST}/triage/upload/foo/bar/test3.txt").times(2)

    expect(Polyphemus::IngestFile.find(name: /test3.txt/)[:ingested_at]).not_to eq(nil)
    expect(Polyphemus::IngestFile.exclude(name: /test3.txt/).map { |f| f[:ingested_at] }).to eq([nil, nil])
  end

  def stub_data
    create(:ingest_file, name: "foo/bar/test1.txt", host: "sftp.example.com", updated_at: "2021-01-01 00:00:00", should_ingest: false)
    create(:ingest_file, name: "foo/bar/test2.txt", host: "sftp.example.com", updated_at: "2015-01-01 00:00:00", should_ingest: false)
    create(:ingest_file, name: "foo/bar/test3.txt", host: "sftp.example.com", updated_at: "1999-01-01 00:00:00", should_ingest: true)
  end

  def stub_ingest_filesystem
    allow_any_instance_of(Polyphemus::SftpIngestMetisTriageFilesEtl).to receive(:ingest_filesystem).and_return(mock_filesystem)
  end

  def mock_filesystem
    mock = Etna::Filesystem::Mock.new do |fname, opts|
      StringIO.new
    end

    mock.with_writeable("foo/bar/test3.txt") do |f|
      f.write(file_contents)
    end

    mock
  end
end
