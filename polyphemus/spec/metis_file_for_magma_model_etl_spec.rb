describe Polyphemus::MetisFileForMagmaModelEtl do
  before(:each) do
    Polyphemus::EtlExecutor.ensure_for_etl(TestMetisFileForMagmaModelEtl)
  end

  let(:etl_command) do
    Polyphemus::EtlCommand.new
  end

  let(:etl_executor) do
    etl_command.subcommands["test_metis_file_for_magma_model_etl"]
  end

  def run_etl_command(*args)
    cmd, args, kwds = etl_command.find_command(*args)
    cmd.execute(*args, **kwds)
  end

  def setup_client(metis_client, magma_client)
    allow(etl_executor.subcommands["run"].etl).to receive(:metis_client).and_return(metis_client)
    allow(etl_executor.subcommands["reset"].etl).to receive(:metis_client).and_return(metis_client)
    allow(etl_executor.subcommands["run"].etl).to receive(:magma_client).and_return(magma_client)
    allow(etl_executor.subcommands["reset"].etl).to receive(:magma_client).and_return(magma_client)
  end

  class TestMetisFileForMagmaModelEtl < Polyphemus::MetisFileForMagmaModelEtl
    def initialize(**args)
      super(
        project_bucket_pairs: [["ipi", "integral_data"]],
        limit: 2,
        model_name: "rna_seq",
        file_name_globs: ["BulkRNASeq/**/*.fastq.gz"],
        **args,
      )
    end

    attr_reader :process_calls

    def process(cursor, files)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, files.dup]
    end
  end

  # To re-record this test, delete the metis_file_for_magma_model_etl.e2e.yml cassette and TOKEN=xxxxx to your environment before running.
  it "should process metis files, and support reset" do
    VCR.use_cassette("metis_file_for_magma_model_etl.e2e") do
      metis_client = Etna::Clients::Metis.new(host: "https://metis.development.local", token: ENV["TOKEN"] || TEST_TOKEN)
      magma_client = Etna::Clients::Magma.new(host: "https://magma.development.local", token: ENV["TOKEN"] || TEST_TOKEN)

      setup_client(metis_client, magma_client)

      etl = etl_executor.subcommands["run"].etl

      run_etl_command("test_metis_file_for_magma_model_etl", "run")
      expect(etl.process_calls.length).to_not eq(0)
      expect(etl.process_calls.length).to eq(1)
      check_process_calls(etl)

      etl.process_calls.clear
      run_etl_command("test_metis_file_for_magma_model_etl", "run")
      expect(etl.process_calls.length).to eq(0)

      run_etl_command("test_metis_file_for_magma_model_etl", "reset")

      etl.process_calls.clear
      run_etl_command("test_metis_file_for_magma_model_etl", "run")
      expect(etl.process_calls.length).to_not eq(0)
      check_process_calls(etl)
    end
  end

  def check_process_calls(etl)
    all_files = []
    etl.process_calls.each do |call|
      files = call.last
      all_files.push(*files)
    end

    # Do not duplicate any yields
    expect(all_files).to eq(all_files.uniq { |f| f.file_paths_hashes })

    # Ensure all data is eventually yielded
    expect(all_files.map(&:file_paths_hashes).flatten).to eq(etl.metis_client.find(Etna::Clients::Metis::FindRequest.new(
      project_name: "ipi",
      bucket_name: "integral_data",
      params: [
        Etna::Clients::Metis::FindParam.new(
          type: "file",
          attribute: "name",
          predicate: "glob",
          value: "BulkRNASeq/**/*.fastq.gz",
        ),
        Etna::Clients::Metis::FindParam.new(
          type: "file",
          attribute: "name",
          predicate: "glob",
          value: "IPIADR001.T1.rna.live/**/*",
        ),
      ],
    )).files.all.map { |f| [f.file_path, f.file_hash] }.flatten)
  end
end
