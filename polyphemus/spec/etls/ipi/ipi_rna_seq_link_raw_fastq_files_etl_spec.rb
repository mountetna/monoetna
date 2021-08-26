describe Polyphemus::IpiRnaSeqLinkRawFastqFilesEtl do
  let(:bucket_name) { "test" }
  let(:project_name) { "ipi" }
  let(:cursor) {
    Polyphemus::MetisFileInWatchFolderCursor.new(
      job_name: "test",
      project_name: project_name,
      bucket_name: bucket_name,
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
      bucket_name: bucket_name,
      project_name: project_name,
    })
  end

  describe "updates Magma records" do
    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files" do
      etl = Polyphemus::IpiRnaSeqLinkRawFastqFilesEtl.new

      etl.process(cursor, [
        file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz"),
        file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz"),
        file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz"),
      ])

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))
    end

    it "when file has been deleted" do
      etl = Polyphemus::IpiRnaSeqLinkRawFastqFilesEtl.new

      etl.process(cursor, [
        file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz"),
        file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz"),
        file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz"),
      ])

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))

      etl.process(cursor, [
        file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz"),
        file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz"),
      ])

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))
    end

    it "when file hash changes" do
      etl = Polyphemus::IpiRnaSeqLinkRawFastqFilesEtl.new

      etl.process(cursor, [
        file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz", "hash1"),
        file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz", "hash2"),
        file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz", "hash3"),
      ])

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))

      etl.process(cursor, [
        file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz", "hash1"),
        file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz", "new-hash2"),
        file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz", "hash3"),
      ])

      # Make sure rna_seq records are updated again
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "PATIENT001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/test/BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "PATIENT001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 })).twice
    end
  end
end
