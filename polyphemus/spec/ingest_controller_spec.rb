describe Polyphemus::Server do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  it "cannot list from invalid host" do
    stub_ingest_files
    auth_header(:administrator)
    get("/test/ingest/list/sftp.foo.com/foo/bar")

    expect(last_response.status).to eq(422)
  end

  it "does not allow non-admins to list directories" do
    stub_ingest_files
    auth_header(:editor)
    get("/test/ingest/list/sftp.foo.com/foo/bar")

    expect(last_response.status).to eq(403)
  end

  it "returns list of files" do
    stub_ingest_files([{
      name: "foo/bar/test1.txt",
      host: "sftp.example.com",
      updated_at: "2021-01-01 00:00:00",
      should_ingest: false,
    }, {
      name: "foo2/bar/test2.txt",
      host: "sftp.example.com",
      updated_at: "2015-01-01 00:00:00",
      should_ingest: false,
    }, {
      name: "foo/bar/test3.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: true,
      ingested_at: "2021-01-01 00:00:00",
    }])

    auth_header(:administrator)
    get("/test/ingest/list/sftp.example.com/foo/bar")

    expect(last_response.status).to eq(200)
    expect(json_body[:results]).to eq([
      {
        name: "foo/bar/test1.txt",
        host: "sftp.example.com",
      },
    ])
  end

  xit "can update should_ingest" do
  end
end
