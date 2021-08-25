describe Polyphemus::MetisFileInWatchFolderEtl do
  before(:each) do
    Polyphemus::EtlExecutor.ensure_for_etl(TestMetisFileInWatchFolderEtl)
  end

  let(:etl_command) do
    Polyphemus::EtlCommand.new
  end

  let(:etl_executor) do
    etl_command.subcommands["test_metis_file_in_watch_folder_etl"]
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

  class TestMetisFileInWatchFolderEtl < Polyphemus::MetisFileInWatchFolderEtl
    def initialize(**args)
      super(
        project_bucket_model_tuples: [["ipi", "integral_data", "rna_seq"]],
        file_name_globs: ["BulkRNASeq/**/*.fastq.gz"],
        metis_path_to_record_name_regex: Polyphemus::IpiRnaSeqLinkFastQsEtl::PATH_REGEX,
        **args,
      )
    end

    attr_reader :process_calls

    def process(cursor, files)
      @process_calls ||= []
      @process_calls << [cursor.updated_at.dup, cursor.value.dup, files.dup]
    end
  end

  it "should process metis files in watch folders, and support reset" do
    stub_watch_folders
    stub_metis_setup
    stub_list_folder(bucket: RELEASE_BUCKET, response_body: {
                       files: [{
                         project_name: PROJECT,
                         bucket_name: RELEASE_BUCKET,
                         file_path: "path1/path1_1/sample.txt",
                       }, {
                         project_name: PROJECT,
                         bucket_name: RELEASE_BUCKET,
                         file_path: "path1/path1_1/path1_1_1/sample.txt",
                       }],
                     }.to_json)

    etl = etl_executor.subcommands["run"].etl

    run_etl_command("test_metis_file_in_watch_folder_etl", "run")
    expect(etl.process_calls.length).to eq(1)

    etl.process_calls.clear
    run_etl_command("test_metis_file_in_watch_folder_etl", "run")
    expect(etl.process_calls.length).to eq(0)

    run_etl_command("test_metis_file_in_watch_folder_etl", "reset")

    etl.process_calls.clear
    run_etl_command("test_metis_file_in_watch_folder_etl", "run")
    expect(etl.process_calls.length).to eq(1)
  end
end
