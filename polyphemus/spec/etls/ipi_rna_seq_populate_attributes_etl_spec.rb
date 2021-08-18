describe Polyphemus::IpiRnaSeqPopulateAttributesEtl do
  before(:each) do
    stub_metis_setup
    @all_updates = []
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
    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files" do
      stub_metis_scan("Plate1", [
        file("rnaseq_table.tsv", "plate1_rnaseq_new/results/rnaseq_table.tsv"),
      ])
      stub_download_file({
        project: "ipi",
        file_contents: ::File.read("spec/fixtures/ipi_rna_seq_results.tsv"),
      })

      etl = Polyphemus::IpiRnaSeqPopulateAttributesEtl.new

      etl.run_once

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "well": "2",
                                         "pipeline_version": "v1.0",
                                         "ribosomal_read_count": 12,
                                         "utr_pct": 29.2,
                                         "filtered_mean_length": 90,
                                       },
                                     },
                                   },
                                 }))
    end
  end

  def stub_metis_scan(record_name, change_list)
    allow_any_instance_of(Polyphemus::HashScanBasedEtlScanner).to receive(:find_batch).and_return(
      [Polyphemus::MetisFilesForMagmaRecord.new(record_name, change_list)]
    )
  end
end
