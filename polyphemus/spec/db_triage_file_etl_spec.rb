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
    stub_ingest_files

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
end
