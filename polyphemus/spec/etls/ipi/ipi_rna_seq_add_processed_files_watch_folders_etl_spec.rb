describe Polyphemus::IpiRnaSeqAddProcessedFilesWatchFoldersEtl do
  let(:project_name) { "ipi" }
  let(:bucket_name) { "data" }
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

  def folder(folder_name, folder_path, updated_at = Time.now)
    Etna::Clients::Metis::Folder.new({
      folder_name: folder_name,
      folder_path: folder_path,
      updated_at: updated_at,
    })
  end

  describe "create Polyphemus::WatchFile records" do
    it "for invalid NASH / NAFLD samples" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::IpiRnaSeqAddProcessedFilesWatchFoldersEtl.new

      etl.process(cursor, [
        folder("IPIADR001.NASH1.rna.live", "plate1_rnaseq_new/output/IPIADR001.NASH1.rna.live"),
        folder("IPIADR001.NAFLD1.rna.live", "plate1_rnaseq_new/output/IPIADR001.NAFLD1.rna.live"),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(2)
    end

    it "for incorrectly named samples" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::IpiRnaSeqAddProcessedFilesWatchFoldersEtl.new

      etl.process(cursor, [
        folder("WRONG001.T1.rna.tumor", "plate1_rnaseq_new/output/WRONG001.T1.rna.tumor"),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(1)
    end
  end

  it "links files in found folders" do
    stub_list_folder(
      project: project_name,
      bucket: bucket_name,
      response_data: {
        files: [{
          file_name: "plate1_rnaseq_new/output/PATIENT1.N1.rna.live/something.fastq.gz",
          bucket_name: bucket_name,
          project_name: project_name,
        }],
        folders: [],
      },
    )

    etl = Polyphemus::IpiRnaSeqAddProcessedFilesWatchFoldersEtl.new

    etl.process(cursor, [
      folder("PATIENT1.N1.rna.live", "plate1_rnaseq_new/output/PATIENT1.N1.rna.live"),
      folder("PATIENT1.T1.rna.live", "plate1_rnaseq_new/output/PATIENT1.T1.rna.live"),
      folder("PATIENT2.T1.rna.live", "plate2_rnaseq_new/output/PATIENT2.T1.rna.live"),
    ])

    # Make sure rna_seq records are updated. Once per folder.
    expect(WebMock).to have_requested(:post, /#{MAGMA_HOST}\/update/).times(3)
  end
end
