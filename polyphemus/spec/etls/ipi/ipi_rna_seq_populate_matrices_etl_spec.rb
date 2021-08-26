describe Polyphemus::IpiRnaSeqPopulateMatricesEtl do
  let(:helper) { IpiHelper.new("lib/etls/renaming/projects/test_renames.json") }

  before(:each) do
    allow(IpiHelper).to receive(:new).and_return(helper)
    stub_metis_setup
    @all_updates = []
    copy_renaming_project
  end

  def file(file_name, file_path, file_hash = SecureRandom.hex, updated_at = Time.now)
    Etna::Clients::Metis::File.new({
      file_name: file_name,
      file_path: file_path,
      updated_at: updated_at,
      file_hash: file_hash,
      project_name: "ipi",
    })
  end

  describe "updates Magma records" do
    let(:cursor) {
      Polyphemus::MetisFileEtlCursor.new(
        job_name: "test",
        project_name: "ipi",
        bucket_name: "test",
      )
    }

    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files" do
      stub_download_file({
        project: "ipi",
        file_contents: ::File.read("spec/fixtures/ipi_gene_counts.tsv"),
      })

      etl = Polyphemus::IpiRnaSeqPopulateMatricesEtl.new

      etl.process(cursor, [
        file("gene_counts_table.tsv", "plate1_rnaseq_new/results/gene_counts_table.tsv"),
      ])

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "gene_counts": [1.2, 0],
                                       },
                                       "RIGHT001.T1.rna.tumor": {
                                         "gene_counts": [2.5, 0],
                                       },
                                     },
                                   },
                                 }))
    end
  end
end
