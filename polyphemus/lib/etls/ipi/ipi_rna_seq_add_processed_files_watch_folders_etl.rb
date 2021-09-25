require_relative "../add_watch_folder_base_etl"
require_relative "./linkers/ipi_rna_seq_processed_files_linker"

class Polyphemus::IpiRnaSeqAddProcessedFilesWatchFoldersEtl < Polyphemus::AddWatchFolderBaseEtl
  PROJECT = "ipi"
  BUCKET = "data"

  def initialize
    super(
      linker: Polyphemus::IpiRnaSeqProcessedFilesLinker.new(project_name: PROJECT, bucket_name: BUCKET),
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_path_regexes: {
        "#{PROJECT}_#{BUCKET}": /^bulkRNASeq\/.*\/output\/.*$/,
      },
      model_name: "rna_seq",
      watch_type: "link_files",
    )
  end
end
