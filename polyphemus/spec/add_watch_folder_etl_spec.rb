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
    it "for found folders that match the regex" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
        watch_type: "link_files",
        folder_path_regexes: {
          "#{PROJECT}_#{RESTRICT_BUCKET}": /plate1/,
        },
      )

      etl.process(cursor, [
        create_metis_folder("PATIENT1.N1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.N1.rna.live", id: 1),
        create_metis_folder("PATIENT1.T1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.T1.rna.live", id: 2),
        create_metis_folder("PATIENT2.T1.rna.live", "bucket/plate2_rnaseq_new/output/PATIENT2.T1.rna.live", id: 3),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(2)
    end

    it "throws exception if no regex provided" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
        watch_type: "link_files",
      )

      expect {
        etl.process(cursor, [
          create_metis_folder("PATIENT1.N1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.N1.rna.live", id: 1),
          create_metis_folder("PATIENT1.T1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.T1.rna.live", id: 2),
          create_metis_folder("PATIENT2.T1.rna.live", "bucket/plate2_rnaseq_new/output/PATIENT2.T1.rna.live", id: 3),
        ])
      }.to raise_error(Polyphemus::EtlError)

      expect(Polyphemus::WatchFolder.count).to eq(0)
    end

    it "can have duplicate folders with different watch types" do
      expect(Polyphemus::WatchFolder.count).to eq(0)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
        watch_type: "link_files",
        folder_path_regexes: {
          "#{PROJECT}_#{RESTRICT_BUCKET}": /plate1/,
        },
      )

      etl.process(cursor, [
        create_metis_folder("PATIENT1.N1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.N1.rna.live", id: 1),
        create_metis_folder("PATIENT1.T1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.T1.rna.live", id: 2),
        create_metis_folder("PATIENT2.T1.rna.live", "bucket/plate2_rnaseq_new/output/PATIENT2.T1.rna.live", id: 3),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(2)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
        watch_type: "process_files",
        folder_path_regexes: {
          "#{PROJECT}_#{RESTRICT_BUCKET}": /plate1/,
        },
      )

      etl.process(cursor, [
        create_metis_folder("PATIENT1.N1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.N1.rna.live", id: 1),
        create_metis_folder("PATIENT1.T1.rna.live", "bucket/plate1_rnaseq_new/output/PATIENT1.T1.rna.live", id: 2),
        create_metis_folder("PATIENT2.T1.rna.live", "bucket/plate2_rnaseq_new/output/PATIENT2.T1.rna.live", id: 3),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(4)
    end

    it "removes folders who no longer match the globs" do
      changed_folder = create(
        :watch_folder,
        project_name: PROJECT,
        bucket_name: RESTRICT_BUCKET,
        updated_at: "2021-01-01 00:00:00",
        folder_path: "path1/path1_1/path1_1_1",
        metis_id: 1,
        watch_type: "link_files",
      )

      expect(Polyphemus::WatchFolder.count).to eq(1)

      etl = Polyphemus::AddWatchFolderBaseEtl.new(
        project_bucket_pairs: [[PROJECT, RESTRICT_BUCKET]],
        model_name: "foo",
        folder_path_regexes: {
          "#{PROJECT}_#{RESTRICT_BUCKET}": /\/path1_1\//,
        },
        watch_type: "link_files",
      )

      etl.process(cursor, [
        create_metis_folder("PATIENT1.N1.rna.live", "path1/path1_2/path2_1_1", id: 1),
        create_metis_folder("PATIENT1.N1.rna.live", "path1/path1_2/path2_1_2", id: 2),
      ])

      expect(Polyphemus::WatchFolder.count).to eq(0)
    end
  end
end
