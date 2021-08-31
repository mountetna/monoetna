describe Polyphemus::IpiRnaSeqLinkProcessedFilesEtl do
  let(:bucket_name) { "test" }
  let(:project_name) { "ipi" }
  let(:cursor) {
    Polyphemus::MetisFileEtlCursor.new(
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
    @all_updates = []
    stub_watch_folders([{
      project_name: project_name,
      bucket_name: bucket_name,
      metis_id: 1,
      folder_path: "bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.comp",
      watch_type: "link_files",
    }])
  end

  describe "updates Magma records" do
    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files for file attribute" do
      etl = Polyphemus::IpiRnaSeqLinkProcessedFilesEtl.new

      etl.process(cursor, [
        create_metis_file("PATIENT001.T1.comp.deduplicated.cram", "bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram"),
      ])

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "cram": {
                                           "path": "metis://ipi/test/bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram",
                                           "original_filename": "PATIENT001.T1.comp.deduplicated.cram",
                                         },
                                       },
                                     },
                                   },
                                 }))
    end

    it "when scanner finds new files for file_collection attribute" do
      etl = Polyphemus::IpiRnaSeqLinkProcessedFilesEtl.new

      etl.process(cursor, [
        create_metis_file("PATIENT001.T1.comp.unmapped.1.fastq.gz", "bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.1.fastq.gz"),
        create_metis_file("PATIENT001.T1.comp.unmapped.2.fastq.gz", "bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.2.fastq.gz"),
      ])

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "unmapped_fastqs": [{
                                           "path": "metis://ipi/test/bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.1.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.unmapped.1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.2.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.unmapped.2.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))
    end

    it "correctly ignores non-cancer files" do
      etl = Polyphemus::IpiRnaSeqLinkProcessedFilesEtl.new

      etl.process(cursor, [
        create_metis_file("PATIENT001.T1.NAFLD.blahblah3.junction", "bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.NAFLD/PATIENT001.T1.NAFLD.blahblah3.junction"),
      ])

      # Make sure rna_seq records are NOT updated
      expect(WebMock).not_to have_requested(:post, /#{MAGMA_HOST}\/update/)
                               .with(body: hash_including({
                                       "revisions": {
                                         "rna_seq": {
                                           "PATIENT001.T1.NAFLD": {
                                             "junction": {
                                               "path": "metis://ipi/test/bulkRNASeq/plate1_blahblah/output/PATIENT001.T1.NAFLD/PATIENT001.T1.NAFLD.blahblah3.junction",
                                               "original_filename": "PATIENT001.T1.NAFLD.blahblah3.junction",
                                             },
                                           },
                                         },
                                       },
                                     }))
    end

    it "correctly renames renamed tube_names" do
      stub_watch_folders([{
        project_name: project_name,
        bucket_name: bucket_name,
        metis_id: 2,
        folder_path: "bulkRNASeq/plate1_blahblah/output/WRONG001.T1.rna.tumor",
        watch_type: "link_files",
      }])

      etl = Polyphemus::IpiRnaSeqLinkProcessedFilesEtl.new

      etl.process(cursor, [
        create_metis_file("WRONG001.T1.rna.tumor.deduplicated.cram.crai", "bulkRNASeq/plate1_blahblah/output/WRONG001.T1.rna.tumor/WRONG001.T1.rna.tumor.deduplicated.cram.crai", folder_id: 2),
      ])

      # Make sure rna_seq records are updated for renamed patient, but pointing to the "wrong" file locations
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "RIGHT001.T1.rna.tumor": {
                                         "cram_index": {
                                           "path": "metis://ipi/test/bulkRNASeq/plate1_blahblah/output/WRONG001.T1.rna.tumor/WRONG001.T1.rna.tumor.deduplicated.cram.crai",
                                           "original_filename": "WRONG001.T1.rna.tumor.deduplicated.cram.crai",
                                         },
                                       },
                                     },
                                   },
                                 }))
    end
  end
end
