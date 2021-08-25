require_relative "../add_watch_folder_base_etl"
require_relative "../../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqPopulateWatchFoldersEtl < Polyphemus::AddWatchFolderBaseEtl
  PROJECT = "ipi"
  BUCKET = "data"
  MODEL = "rna_seq"

  def initialize
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_name_globs: ["output/*", "bulkRNASeq/*"],
    )
  end

  def should_skip_folder?(folder)
    !current_magma_record_names.include?(folder.folder_name)
  end

  def current_magma_record_names
    @current_magma_record_names ||= begin
        request = Etna::Clients::Magma::RetrievalRequest.new(
          project_name: PROJECT,
          model_name: MODEL,
          attribute_names: ["identifier"],
          record_names: "all",
          hide_templates: true,
        )
        self.magma_client.retrieve(request).models.model(MODEL).documents.document_keys
      end
  end
end
