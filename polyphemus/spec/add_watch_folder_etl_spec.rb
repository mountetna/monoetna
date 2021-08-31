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

  describe "create Polyphemus::WatchFile records" do
    it "for found folders" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
        watch_type: "link_files",
      )

      etl.process(cursor, [
        create_metis_folder("PATIENT1.N1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.N1.rna.live", id: 1),
        create_metis_folder("PATIENT1.T1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.T1.rna.live", id: 2),
        create_metis_folder("PATIENT2.T1.rna.live", "bucket/plate2_rnaseq_new/output/PATIENT2.T1.rna.live", id: 3),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(3)
    end

    it "removes folders who no longer match the globs" do
      changed_folder = create(
        :watch_folder,
        project_name: PROJECT,
        bucket_name: RESTRICT_BUCKET,
        updated_at: "2021-01-01 00:00:00",
        folder_path: "path1/path1_1/path1_1_1",
        folder_id: 1,
        watch_type: "link_files",
      )

      expect(Polyphemus::WatchFolder.count).to eq(1)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
        folder_name_globs: ["path1_1/*", "path1/*"],
        watch_type: "link_files",
      )

      etl.process(cursor, [
        create_metis_folder("PATIENT1.N1.rna.live", "path1/path1_2/path1_1_1", id: 1),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(0)
    end
  end
end
