describe Polyphemus::SyncCatFilesEtl do
  before(:each) do
    stub_data
  end

  it "should add new, uningested records that match our oligo" do
    stub_rsync_data([
      MockChange.new("foo/bar/oligo-test1.txt"),
      MockChange.new("foo/bar/oligo-test2.txt"),
      MockChange.new("foo/bar/oligo-test3.txt"),
      MockChange.new("wrong-olino-file.txt"),
      MockChange.new("foo/oligo/file-to-add.txt"),
    ])
    sync_etl = Polyphemus::SyncCatFilesEtl.new

    expect(Polyphemus::IngestFile.count).to eq(3)

    sync_etl.run_once

    expect(Polyphemus::IngestFile.count).to eq(4)
    expect(Polyphemus::IngestFile.last[:should_ingest]).to eq(true)
    expect(Polyphemus::IngestFile.last[:name]).to eq("foo/oligo/file-to-add.txt")
  end

  it "should remove old, uningested records" do
    stub_rsync_data([
      MockChange.new("foo/bar/oligo-test1.txt"),
      MockChange.new("foo/bar/oligo-test3.txt"),
    ])
    sync_etl = Polyphemus::SyncCatFilesEtl.new

    expect(Polyphemus::IngestFile.count).to eq(3)
    expect(Polyphemus::IngestFile.where(removed_from_source: false).count).to eq(3)

    sync_etl.run_once

    expect(Polyphemus::IngestFile.count).to eq(3)
    expect(Polyphemus::IngestFile.where(removed_from_source: false).count).to eq(2)
  end

  it "should not change ingested records" do
    stub_rsync_data([
      MockChange.new("foo/bar/oligo-test1.txt"),
      MockChange.new("foo/bar/oligo-test2.txt"),
      MockChange.new("foo/bar/oligo-test3.txt"),
    ])
    sync_etl = Polyphemus::SyncCatFilesEtl.new

    expect(Polyphemus::IngestFile.count).to eq(3)
    test3 = Polyphemus::IngestFile.find(name: /test3.txt/)

    expect(test3[:triage_ingested_at]).to eq(time_at("2021-01-01 00:00:00"))
    expect(test3[:should_ingest]).to eq(true)

    sync_etl.run_once

    expect(Polyphemus::IngestFile.count).to eq(3)
    expect(test3[:triage_ingested_at]).to eq(time_at("2021-01-01 00:00:00"))
    expect(test3[:should_ingest]).to eq(true)
  end

  def stub_data
    stub_ingest_files([{
      name: "foo/bar/oligo-test1.txt",
      host: "sftp.example.com",
      updated_at: "2021-01-01 00:00:00",
      should_ingest: false,
    }, {
      name: "foo/bar/oligo-test2.txt",
      host: "sftp.example.com",
      updated_at: "2015-01-01 00:00:00",
      should_ingest: false,
    }, {
      name: "foo/bar/oligo-test3.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: true,
      triage_ingested_at: "2021-01-01 00:00:00",
    }])
  end

  def stub_rsync_data(change_list)
    allow_any_instance_of(Polyphemus::FilenameScanBasedEtlScanner).to receive(:find_batch).and_return(change_list)
  end

  def time_at(time_string)
    Time.at(DateTime.parse(time_string).strftime("%s").to_i)
  end

  class MockChange
    def initialize(filename)
      @filename = filename
    end

    def filename
      @filename
    end

    def timestamp
      "new"
    end
  end
end
