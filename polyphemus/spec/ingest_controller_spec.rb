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

  it "returns list of uningested / un-queued files, including subdirectories" do
    stub_ingest_files([{
      name: "foo/bar/test1.txt",
      host: "sftp.example.com",
      updated_at: "2021-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: false,
    }, {
      name: "foo2/bar/test2.txt",
      host: "sftp.example.com",
      updated_at: "2015-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: false,
    }, {
      name: "foo/bar/test3.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: true,
      ingested_at: "2021-01-01 00:00:00",
      removed_from_source: false,
    }, {
      name: "foo/bar2/test4.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: false,
    }, {
      name: "foo/bar/test3/test5.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: false,
    }, {
      name: "foo/bar/test3/test6.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: true,
    }])

    auth_header(:administrator)
    get("/test/ingest/list/sftp.example.com/foo/bar")

    expect(last_response.status).to eq(200)
    expect(json_body[:files]).to eq([
      {
        name: "foo/bar/test1.txt",
        host: "sftp.example.com",
      }, {
        name: "foo/bar/test3/test5.txt",
        host: "sftp.example.com",
      },
    ])

    get("/test/ingest/list/sftp.example.com/foo/bar/")

    expect(last_response.status).to eq(200)
    expect(json_body[:files]).to eq([
      {
        name: "foo/bar/test1.txt",
        host: "sftp.example.com",
      }, {
        name: "foo/bar/test3/test5.txt",
        host: "sftp.example.com",
      },
    ])
  end

  it "can queue up files for ingestion" do
    stub_ingest_files([{
      name: "foo/bar/test1.txt",
      host: "sftp.example.com",
      updated_at: "2021-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: false,
    }, {
      name: "foo2/bar/test2.txt",
      host: "sftp.example.com",
      updated_at: "2015-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: false,
    }, {
      name: "foo/bar/test3.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: true,
      ingested_at: "2021-01-01 00:00:00",
      removed_from_source: false,
    }, {
      name: "foo/bar2/test4.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: false,
      removed_from_source: false,
    }, {
      name: "foo/bar/test3/test5.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: true,
      removed_from_source: false,
    }])

    auth_header(:administrator)
    get("/test/ingest/list/sftp.example.com/foo/bar")

    expect(last_response.status).to eq(200)
    expect(json_body[:files]).to eq([{
                                   name: "foo/bar/test1.txt",
                                   host: "sftp.example.com",
                                 }])

    post("/test/ingest/enqueue/sftp.example.com/foo/bar")

    expect(last_response.status).to eq(200)

    get("/test/ingest/list/sftp.example.com/foo/bar")

    expect(last_response.status).to eq(200)
    expect(json_body[:files]).to eq([])
  end
end
