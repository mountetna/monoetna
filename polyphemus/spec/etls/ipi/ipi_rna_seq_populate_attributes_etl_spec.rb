describe Polyphemus::IpiRnaSeqPopulateAttributesEtl do
  let(:bucket_name) { "test" }
  let(:project_name) { "ipi" }

  before(:each) do
    stub_metis_setup
    @all_updates = []
    stub_watch_folders([{
      project_name: project_name,
      bucket_name: bucket_name,
      folder_id: 1,
      folder_path: "plate1_rnaseq_new/results/",
      watch_type: "process_files",
    }])
  end

  describe "updates Magma records" do
    let(:cursor) {
      Polyphemus::MetisFileEtlCursor.new(
        job_name: "test",
        project_name: project_name,
        bucket_name: bucket_name,
      )
    }

    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files" do
      stub_download_file({
        project: project_name,
        file_contents: ::File.read("spec/fixtures/ipi_rna_seq_results.tsv"),
      })

      etl = Polyphemus::IpiRnaSeqPopulateAttributesEtl.new

      etl.process(cursor, [
        create_metis_file("rnaseq_table.tsv", "plate1_rnaseq_new/results/rnaseq_table.tsv"),
      ])

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
