describe Polyphemus::DbTriageFileNotificationEtl do
  before(:each) do
    Polyphemus::EtlExecutor.ensure_for_etl(TestDbTriageFileNotificationEtl)
  end

  let(:etl_command) do
    Polyphemus::EtlCommand.new
  end

  let(:etl_executor) do
    etl_command.subcommands["test_db_triage_file_notification_etl"]
  end

  def run_etl_command(*args)
    cmd, args, kwds = etl_command.find_command(*args)
    cmd.execute(*args, **kwds)
  end

  class TestDbTriageFileNotificationEtl < Polyphemus::DbTriageFileNotificationEtl
    def initialize(**args)
      super(
        project_bucket_pairs: [["triage", "waiting_room"]],
        limit: 2,
        column_name: :triage_ingested_at,
        **args,
      )
    end

    attr_reader :process_calls

    def process(cursor, files)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, files.dup]
    end
  end

  it "should process c4 ingested files from the given triage table" do
    ingested_files = stub_ingest_files([{
      name: "foo/bar/oligo-test1.txt",
      host: "sftp.example.com",
      updated_at: "2021-01-01 00:00:00",
      should_ingest: false,
      archive_ingested_at: "2021-01-01 00:00:00",
    }, {
      name: "foo/bar/oligo-test2.txt",
      host: "sftp.example.com",
      updated_at: "2015-01-01 00:00:00",
      should_ingest: true,
      triage_ingested_at: "2021-01-01 00:00:00",
    }, {
      name: "foo/bar/oligo-test3.txt",
      host: "sftp.example.com",
      updated_at: "1999-01-01 00:00:00",
      should_ingest: true,
      triage_ingested_at: "2021-01-01 00:00:00",
    }])

    etl = etl_executor.subcommands["run"].etl

    run_etl_command("test_db_triage_file_notification_etl", "run")
    expect(etl.process_calls.length).to eq(1)
    expect(etl.process_calls.first.last.length).to eq(2)

    etl.process_calls.clear
    run_etl_command("test_db_triage_file_notification_etl", "run")
    expect(etl.process_calls.length).to eq(0)

    run_etl_command("test_db_triage_file_notification_etl", "reset")

    etl.process_calls.clear
    run_etl_command("test_db_triage_file_notification_etl", "run")
    expect(etl.process_calls.length).to eq(1)
    expect(etl.process_calls.first.last.length).to eq(2)

    etl.process_calls.clear
    run_etl_command("test_db_triage_file_notification_etl", "run")
    expect(etl.process_calls.length).to eq(0)

    # simulate 'removing' files from the cat and updating the state.
    ingested_files.each do |f|
      # Ensure the update is in the future to prevent race conditions in the test.
      f.update(removed_from_source: true, updated_at: Time.now + 10.seconds)
    end

    etl.process_calls.clear
    run_etl_command("test_db_triage_file_notification_etl", "run")
    expect(etl.process_calls.length).to eq(0)
  end
end
