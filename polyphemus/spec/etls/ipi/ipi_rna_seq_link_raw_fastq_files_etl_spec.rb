describe Polyphemus::IpiRnaSeqLinkRawFastqFilesEtl do
  let(:bucket_name) { "test" }
  let(:project_name) { "ipi" }
  let(:cursor) {
    Polyphemus::MetisFileEtlCursor.new(
      job_name: "test",
      project_name: project_name,
      bucket_name: bucket_name,
    )
  }

  before(:each) do
    stub_metis_setup
    @all_updates = []
    stub_watch_folders([{
      project_name: project_name,
      bucket_name: bucket_name,
      metis_id: 1,
      folder_path: "BulkRNASeq/PATIENT001.T1.comp",
      watch_type: "link_files",
    }])
  end

  describe "updates Magma records" do
    before(:each) do
      stub_magma_update_json
      stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
    end

    it "when scanner finds new files" do
      etl = Polyphemus::IpiRnaSeqLinkRawFastqFilesEtl.new

      etl.process(cursor, [
        create_metis_file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz"),
        create_metis_file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz"),
        create_metis_file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz"),
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
        create_metis_file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz"),
        create_metis_file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz"),
        create_metis_file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz"),
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
        create_metis_file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz"),
        create_metis_file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz"),
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
        create_metis_file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz", file_hash: "hash1"),
        create_metis_file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz", file_hash: "hash2"),
        create_metis_file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz", file_hash: "hash3"),
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
        create_metis_file("PATIENT001.T1.comp.blahblah1.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah1.fastq.gz", file_hash: "hash1"),
        create_metis_file("PATIENT001.T1.comp.blahblah2.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah2.fastq.gz", file_hash: "new-hash2"),
        create_metis_file("PATIENT001.T1.comp.blahblah3.fastq.gz", "BulkRNASeq/PATIENT001.T1.comp/PATIENT001.T1.comp.blahblah3.fastq.gz", file_hash: "hash3"),
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
