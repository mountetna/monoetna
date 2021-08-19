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
    let(:cursor) {
      Polyphemus::MetisFileForMagmaModelEtlCursor.new(
        job_name: "test",
        project_name: "ipi",
        bucket_name: "test",
        model_name: "test",
      )
    }

    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files" do
      stub_download_file({
        project: "ipi",
        file_contents: ::File.read("spec/fixtures/ipi_rna_seq_results.tsv"),
      })

      etl = Polyphemus::IpiRnaSeqPopulateAttributesEtl.new

      etl.process(cursor, mock_metis_files_for_record_scan(
        "Plate1", [
          file("rnaseq_table.tsv", "plate1_rnaseq_new/results/rnaseq_table.tsv"),
        ]
      ))

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
                                         "compartment": "comp2",
                                       },
                                       "PATIENT002.T1.comp": {
                                         "well": "3",
                                         "pipeline_version": "v1.1",
                                         "ribosomal_read_count": 13,
                                         "utr_pct": 99.9,
                                         "filtered_mean_length": 30,
                                         "compartment": "other",
                                       },
                                       "PATIENT003.T1.comp": {
                                         "well": "4",
                                         "pipeline_version": "v1.7",
                                         "ribosomal_read_count": 14,
                                         "utr_pct": 50.0,
                                         "filtered_mean_length": 10,
                                         "compartment": "comp1",
                                       },
                                     },
                                   },
                                 }))
    end
  end
end
