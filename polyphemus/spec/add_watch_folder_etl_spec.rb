describe Polyphemus::AddWatchFolderBaseEtl do
  let(:cursor) {
    Polyphemus::MetisFolderEtlCursor.new(
      job_name: "test",
      project_name: PROJECT,
      bucket_name: RESTRICT_BUCKET,
    )
  }

  before(:each) do
    stub_metis_setup
    stub_list_folder
  end

  def folder(folder_name, folder_path, updated_at = Time.now)
    Etna::Clients::Metis::Folder.new({
      folder_name: folder_name,
      folder_path: folder_path,
      updated_at: updated_at,
    })
  end

  describe "create Polyphemus::WatchFile records" do
    it "for found folders" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
      )

      etl.process(cursor, [
        folder("PATIENT1.N1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.N1.rna.live"),
        folder("PATIENT1.T1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.T1.rna.live"),
        folder("PATIENT2.T1.rna.live", "bucket/plate2_rnaseq_new/output/PATIENT2.T1.rna.live"),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(3)
    end
  end
end
