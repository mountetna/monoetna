describe Polyphemus::MetisPropagateFolderUpdatedAtEtl do
  let(:project_name) { "ipi" }
  let(:bucket_name) { "data" }
  let(:cursor) {
    Polyphemus::MetisFolderEtlCursor.new(
      job_name: "test",
      project_name: project_name,
      bucket_name: bucket_name,
    )
  }

  before(:each) do
    stub_metis_setup
  end

  describe "Polyphemus::MetisPropagateFolderUpdatedAtEtl" do
    it "touches children folders with updated_at older than parent" do
      now = DateTime.now
      Timecop.freeze(now - 100)

      stub_touch_folder(
        project: project_name,
        bucket: bucket_name,
      )
      stub_list_folder(
        url_verb: "list_by_id",
        project: project_name,
        bucket: bucket_name,
        response_body: {
          files: [],
          folders: [{
            project_name: project_name,
            bucket_name: bucket_name,
            folder_path: "output/PATIENT001.T1.comp",
            updated_at: (now - 200).iso8601,
          }],
        },
      )

      etl = Polyphemus::MetisPropagateFolderUpdatedAtEtl.new

      etl.process(cursor, [
        create_metis_folder("output", "output", updated_at: now - 100, id: 1, project_name: project_name, bucket_name: bucket_name),
      ])

      expect(WebMock).to have_requested(:get, /#{METIS_HOST}\/#{project_name}\/folder\/touch\/#{bucket_name}\/output\/PATIENT001.T1.comp/)

      Timecop.return
    end

    it "does not touch children folders with updated_at newer than parent" do
      now = DateTime.now
      Timecop.freeze(now - 100)

      stub_list_folder(
        url_verb: "list_by_id",
        project: project_name,
        bucket: bucket_name,
        response_body: {
          files: [],
          folders: [{
            project_name: project_name,
            bucket_name: bucket_name,
            folder_path: "output/PATIENT001.T1.comp",
            updated_at: (now - 10).iso8601,
          }],
        },
      )

      etl = Polyphemus::MetisPropagateFolderUpdatedAtEtl.new

      etl.process(cursor, [
        create_metis_folder("output", "output", updated_at: now - 100, id: 1, project_name: project_name, bucket_name: bucket_name),
      ])

      expect(WebMock).not_to have_requested(:get, /#{METIS_HOST}\/#{project_name}\/folder\/touch\/#{bucket_name}\/output\/PATIENT001.T1.comp/)

      Timecop.return
    end
  end
end
