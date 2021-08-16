describe Polyphemus::IpiRnaSeqLinkFastQsEtl do
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
      stub_metis_scan("IPIADR001.T1.comp", [
        file("IPIADR001.T1.comp.blahblah1.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz"),
        file("IPIADR001.T1.comp.blahblah2.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz"),
        file("IPIADR001.T1.comp.blahblah3.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz"),
      ])

      etl = Polyphemus::IpiRnaSeqLinkFastQsEtl.new

      etl.run_once

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "IPIADR001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))
    end

    it "when file has been deleted" do
      stub_metis_scan("IPIADR001.T1.comp", [
        file("IPIADR001.T1.comp.blahblah1.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz"),
        file("IPIADR001.T1.comp.blahblah2.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz"),
        file("IPIADR001.T1.comp.blahblah3.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz"),
      ])

      etl = Polyphemus::IpiRnaSeqLinkFastQsEtl.new

      etl.run_once

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "IPIADR001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))

      stub_metis_scan("IPIADR001.T1.comp", [
        file("IPIADR001.T1.comp.blahblah1.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz"),
        file("IPIADR001.T1.comp.blahblah3.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz"),
      ])

      etl.run_once

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "IPIADR001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))
    end

    it "when file hash changes" do
      stub_metis_scan("IPIADR001.T1.comp", [
        file("IPIADR001.T1.comp.blahblah1.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz", "hash1"),
        file("IPIADR001.T1.comp.blahblah2.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz", "hash2"),
        file("IPIADR001.T1.comp.blahblah3.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz", "hash3"),
      ])

      etl = Polyphemus::IpiRnaSeqLinkFastQsEtl.new

      etl.run_once

      # Make sure rna_seq records are updated
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "IPIADR001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 }))

      stub_metis_scan("IPIADR001.T1.comp", [
        file("IPIADR001.T1.comp.blahblah1.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz", "hash1"),
        file("IPIADR001.T1.comp.blahblah2.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz", "new-hash2"),
        file("IPIADR001.T1.comp.blahblah3.fastq.gz", "IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz", "hash3"),
      ])

      etl.run_once

      # Make sure rna_seq records are again
      expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
                           .with(body: hash_including({
                                   "revisions": {
                                     "rna_seq": {
                                       "IPIADR001.T1.comp": {
                                         "raw_fastq_files": [{
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah1.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah1.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah2.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah2.fastq.gz",
                                         }, {
                                           "path": "metis://ipi/integral_data/IPIADR001.T1.comp/IPIADR001.T1.comp.blahblah3.fastq.gz",
                                           "original_filename": "IPIADR001.T1.comp.blahblah3.fastq.gz",
                                         }],
                                       },
                                     },
                                   },
                                 })).twice
    end
  end

  def stub_metis_scan(record_name, change_list)
    allow_any_instance_of(Polyphemus::HashScanBasedEtlScanner).to receive(:find_batch).and_return(
      [Polyphemus::MetisFilesForMagmaRecord.new(record_name, change_list)]
    )
  end
end
