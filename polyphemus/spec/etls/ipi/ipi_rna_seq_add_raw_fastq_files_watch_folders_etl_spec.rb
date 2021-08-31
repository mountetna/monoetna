describe Polyphemus::IpiRnaSeqAddRawFastqFilesWatchFoldersEtl do
  let(:project_name) { "ipi" }
  let(:bucket_name) { "integral_data" }
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
    stub_list_folder(project: project_name, bucket: bucket_name)
    stub_magma_update_json
    stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
  end

  describe "create Polyphemus::WatchFile records" do
    it "for invalid NASH / NAFLD samples" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::IpiRnaSeqAddRawFastqFilesWatchFoldersEtl.new

      etl.process(cursor, [
        create_metis_folder("IPIADR001.NASH1.rna.live", "some_folder/BulkRNASeq/IPIADR001.NASH1.rna.live", id: 1),
        create_metis_folder("IPIADR001.NAFLD1.rna.live", "some_folder/BulkRNASeq/IPIADR001.NAFLD1.rna.live", id: 2),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(2)
    end

    it "for incorrectly named samples" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::IpiRnaSeqAddRawFastqFilesWatchFoldersEtl.new

      etl.process(cursor, [
        create_metis_folder("WRONG001.T1.rna.tumor", "some_folder/BulkRNASeq/WRONG001.T1.rna.tumor", id: 1),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(1)
    end
  end

  it "links files in found folders" do
    stub_list_folder(
      project: project_name,
      bucket: bucket_name,
      response_body: {
        files: [{
          file_name: "something.fastq.gz",
          folder_id: 1,
          bucket_name: bucket_name,
          project_name: project_name,
        }],
        folders: [],
      },
    )

    etl = Polyphemus::IpiRnaSeqAddRawFastqFilesWatchFoldersEtl.new

    etl.process(cursor, [
      create_metis_folder("PATIENT001.T1.comp", "BulkRNASeq/PATIENT001.T1.comp", id: 1),
      create_metis_folder("PATIENT001.N1.comp", "BulkRNASeq/PATIENT001.N1.comp", id: 2),
      create_metis_folder("PATIENT002.T1.comp", "BulkRNASeq/PATIENT002.T1.comp", id: 3),
    ])

    # Make sure rna_seq records are updated. Once per folder with files.
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/)
  end
end
