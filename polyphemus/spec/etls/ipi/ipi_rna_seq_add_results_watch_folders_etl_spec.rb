describe Polyphemus::IpiRnaSeqAddResultsWatchFoldersEtl do
  let(:project_name) { "ipi" }
  let(:bucket_name) { "data" }
  let(:cursor) {
    Polyphemus::MetisFolderEtlCursor.new(
      job_name: "test",
      project_name: project_name,
      bucket_name: bucket_name,
    )
  }
  let(:helper) { IpiHelper.new("lib/etls/renaming/projects/test_renames.json") }

  before(:each) do
    allow(IpiHelper).to receive(:new).and_return(helper)
    stub_metis_setup
    copy_renaming_project
    stub_list_folder(
      url_verb: "list_by_id",
      project: project_name,
      bucket: bucket_name,
    )
    stub_magma_update_json
    stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
  end

  describe "create Polyphemus::WatchFile records" do
    it "for invalid NASH / NAFLD samples" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::IpiRnaSeqAddResultsWatchFoldersEtl.new

      etl.process(cursor, [
        create_metis_folder("results", "bulkRNASeq/plate1_rnaseq_new/results", id: 1),
        create_metis_folder("results", "bulkRNASeq/plate2_rnaseq_new/results", id: 2),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(2)
    end
  end

  describe "processes files in found folders" do
    it "for attributes" do
      stub_download_file(
        project: project_name,
        file_contents: ::File.read("./spec/fixtures/ipi_rna_seq_results.tsv"),
      )
      stub_list_folder(
        url_verb: "list_by_id",
        project: project_name,
        bucket: bucket_name,
        response_body: {
          files: [{
            file_name: "rnaseq_table.tsv",
            download_url: "/#{project_name}/download/rnaseq_table.tsv",
            bucket_name: bucket_name,
            project_name: project_name,
          }],
          folders: [],
        },
      )

      etl = Polyphemus::IpiRnaSeqAddResultsWatchFoldersEtl.new

      etl.process(cursor, [
        create_metis_folder("results", "bulkRNASeq/plate1_rnaseq_new/results", id: 1),
        create_metis_folder("results", "bulkRNASeq/plate2_rnaseq_new/results", id: 2),
      ])

      # Make sure rna_seq records are updated, once per rna_seq record. Only when file is found.
      #   Three valid records in the given response, so three updates.
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).times(3)
    end

    it "for matrices" do
      stub_download_file(
        project: project_name,
        file_contents: ::File.read("./spec/fixtures/ipi_gene_counts.tsv"),
      )
      stub_list_folder(
        url_verb: "list_by_id",
        project: project_name,
        bucket: bucket_name,
        response_body: {
          files: [{
            file_name: "gene_counts_table.tsv",
            download_url: "/#{project_name}/download/gene_counts_table.tsv",
            bucket_name: bucket_name,
            project_name: project_name,
          }],
          folders: [],
        },
      )

      etl = Polyphemus::IpiRnaSeqAddResultsWatchFoldersEtl.new

      etl.process(cursor, [
        create_metis_folder("results", "bulkRNASeq/plate1_rnaseq_new/results", id: 1),
        create_metis_folder("results", "bulkRNASeq/plate2_rnaseq_new/results", id: 2),
      ])

      # Make sure rna_seq records are updated, once per rna_seq record. Only when file is found.
      #   two records in the given response, so two updates.
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).times(2)
    end
  end
end
