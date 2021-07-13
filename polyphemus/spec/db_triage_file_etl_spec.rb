describe Polyphemus::DbTriageFileEtl do
  before(:each) do
    Polyphemus::EtlExecutor.ensure_for_etl(TestDbTriageFileEtl)
  end

  let(:etl_command) do
    Polyphemus::EtlCommand.new
  end

  let(:etl_executor) do
    etl_command.subcommands["test_db_triage_file_etl"]
  end

  def run_etl_command(*args)
    cmd, args, kwds = etl_command.find_command(*args)
    cmd.execute(*args, **kwds)
  end

  class TestDbTriageFileEtl < Polyphemus::DbTriageFileEtl
    def initialize(**args)
      super(
        project_bucket_pairs: [["triage", "waiting_room"]],
        limit: 2,
        table_name: "ingest_files",
        **args,
      )
    end

    attr_reader :process_calls

    def process(cursor, files)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, files.dup]
    end
  end

  it "should process files from the given triage table" do
    stub_data

    etl = etl_executor.subcommands["run"].etl

    run_etl_command("test_db_triage_file_etl", "run")
    expect(etl.process_calls.length).to eq(1)

    etl.process_calls.clear
    run_etl_command("test_db_triage_file_etl", "run")
    expect(etl.process_calls.length).to eq(0)

    run_etl_command("test_db_triage_file_etl", "reset")

    etl.process_calls.clear
    run_etl_command("test_db_triage_file_etl", "run")
    expect(etl.process_calls.length).to eq(1)
  end

  def stub_data
    create(:ingest_file, name: "foo/bar/test1.txt", host: "sftp.example.com", updated_at: "2021-01-01 00:00:00", should_ingest: false)
    create(:ingest_file, name: "foo/bar/test2.txt", host: "sftp.example.com", updated_at: "2015-01-01 00:00:00", should_ingest: false)
    create(:ingest_file, name: "foo/bar/test3.txt", host: "sftp.example.com", updated_at: "1999-01-01 00:00:00", should_ingest: true)
  end
end
