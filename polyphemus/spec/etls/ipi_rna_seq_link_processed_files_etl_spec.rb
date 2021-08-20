describe Polyphemus::IpiRnaSeqLinkProcessedFilesEtl do
  let(:cursor) {
    Polyphemus::MetisFileForMagmaModelEtlCursor.new(
      job_name: "test",
      project_name: "ipi",
      bucket_name: "test",
      model_name: "rna_seq",
    )
  }

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
    })
  end

  describe "updates Magma records" do
    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files" do
      etl = Polyphemus::IpiRnaSeqLinkProcessedFilesEtl.new

      etl.process(cursor, mock_metis_files_for_record_scan(
        "PATIENT001.T1.comp", [
          file("PATIENT001.T1.comp.unmapped.1.fastq.gz", "bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.1.fastq.gz"),
          file("PATIENT001.T1.comp.unmapped.2.fastq.gz", "bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.2.fastq.gz"),
          file("PATIENT001.T1.comp.blahblah3.junction", "bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.junction"),
          file("PATIENT001.T1.comp.deduplicated.cram", "bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram"),
          file("PATIENT001.T1.comp.deduplicated.cram.crai", "bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram.crai"),
        ]
      ))

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "unmapped_fastqs": [{
                                           "path": "metis://ipi/data/bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.1.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.unmapped.1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/data/bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.unmapped.2.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.unmapped.2.fastq.gz",
                                         }],
                                         "cram": {
                                           "path": "metis://ipi/data/bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram",
                                           "original_filename": "PATIENT001.T1.comp.deduplicated.cram",
                                         },
                                         "cram_index": {
                                           "path": "metis://ipi/data/bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram.crai",
                                           "original_filename": "PATIENT001.T1.comp.deduplicated.cram.crai",
                                         },
                                         "junction": {
                                           "path": "metis://ipi/data/bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.junction",
                                           "original_filename": "PATIENT001.T1.comp.blahblah3.junction",
                                         },
                                       },
                                     },
                                   },
                                 }))
    end

    it "sends ::blank or [] when no file found" do
      etl = Polyphemus::IpiRnaSeqLinkProcessedFilesEtl.new

      etl.process(cursor, mock_metis_files_for_record_scan(
        "PATIENT001.T1.comp", [
          file("PATIENT001.T1.comp.deduplicated.cram", "bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram"),
          file("PATIENT001.T1.comp.deduplicated.cram.crai", "bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram.crai"),
        ]
      ))

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "cram": {
                                           "path": "metis://ipi/data/bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram",
                                           "original_filename": "PATIENT001.T1.comp.deduplicated.cram",
                                         },
                                         "cram_index": {
                                           "path": "metis://ipi/data/bulkRNASeq/output/PATIENT001.T1.comp/PATIENT001.T1.comp.deduplicated.cram.crai",
                                           "original_filename": "PATIENT001.T1.comp.deduplicated.cram.crai",
                                         },
                                         "unmapped_fastqs": [],
                                         "junction": {
                                           "path": "::blank",
                                         },
                                       },
                                     },
                                   },
                                 }))
    end
  end
end
