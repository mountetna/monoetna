describe Polyphemus::MagmaRecordFileEtl do
  before(:each) do
    Polyphemus::EtlExecutor.ensure_for_etl(TestMagmaFileEtl)
  end

  let(:etl_command) do
    Polyphemus::EtlCommand.new
  end

  let(:etl_executor) do
    etl_command.subcommands["test_magma_file_etl"]
  end

  def run_etl_command(*args)
    cmd, args, kwds = etl_command.find_command(*args)
    cmd.execute(*args, **kwds)
  end

  def setup_client(magma_client)
    allow(etl_executor.subcommands["run"].etl).to receive(:magma_client).and_return(magma_client)
    allow(etl_executor.subcommands["reset"].etl).to receive(:magma_client).and_return(magma_client)
  end

  class TestMagmaFileEtl < Polyphemus::MagmaRecordFileEtl
    def initialize(**args)
      super(project_model_pairs: [["ipi", "patient"]],
            limit: 2,
            attribute_names: ["updated_at", "flojo_file_processed"],
            **args)
    end

    attr_reader :process_calls

    def process(cursor, records)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, records.dup]
    end
  end

  # To re-record this test, delete the magma_record_file_etl.e2e.yml cassette and TOKEN=xxxxx to your environment before running.
  it "should process magma records, and support reset" do
    VCR.use_cassette("magma_record_file_etl.e2e") do
      magma_client = Etna::Clients::Magma.new(host: "https://magma.development.local", token: ENV["TOKEN"] || TEST_TOKEN)
      setup_client(magma_client)

      etl = etl_executor.subcommands["run"].etl

      run_etl_command("test_magma_file_etl", "run")

      expect(etl.process_calls.length).to_not eq(0)
      expect(etl.process_calls.length).to_not eq(1)
      check_process_calls(etl)

      etl.process_calls.clear
      run_etl_command("test_magma_file_etl", "run")
      expect(etl.process_calls.length).to eq(0)

      run_etl_command("test_magma_file_etl", "reset")

      etl.process_calls.clear
      run_etl_command("test_magma_file_etl", "run")
      expect(etl.process_calls.length).to_not eq(0)
      check_process_calls(etl)
    end
  end

  def check_process_calls(etl)
    all_processed_records = []
    etl.process_calls.each do |call|
      records = call.last
      all_processed_records.push(*records)
    end

    # Do not duplicate any yields
    expect(all_processed_records).to eq(all_processed_records.uniq)

    # Ensure all data is eventually yielded
    remote_records_ids = etl.magma_client.query(Etna::Clients::Magma::QueryRequest.new(
      project_name: "ipi",
      query: ["patient",
              ["::has", "flojo_file_processed"],
              "::all",
              "::identifier"],
    )).answer.map { |r| r.first }.sort

    expect(all_processed_records.map { |r| r.keys.first }.sort).to eq(remote_records_ids)
  end
end
