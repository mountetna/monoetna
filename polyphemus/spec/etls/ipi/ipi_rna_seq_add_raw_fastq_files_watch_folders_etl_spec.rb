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
    stub_list_folder(
      url_verb: "list_by_id",
      project: project_name,
      bucket: bucket_name,
    )
    stub_magma_update_json
    stub_magma_models(fixture: "spec/fixtures/magma_ipi_models_with_records.json")
  end

end
