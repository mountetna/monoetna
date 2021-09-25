require_relative "../add_watch_folder_base_etl"
require_relative "./linkers/ipi_rna_seq_raw_fastq_files_linker"

class Polyphemus::IpiRnaSeqAddRawFastqFilesWatchFoldersEtl < Polyphemus::AddWatchFolderBaseEtl
  PROJECT = "ipi"
  BUCKET = "integral_data"

  def initialize
    super(
      linker: Polyphemus::IpiRnaSeqRawFastqFilesLinker.new(
        project_name: PROJECT, bucket_name: BUCKET,
      ),
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_path_regexes: {
        "#{PROJECT}_#{BUCKET}": /^.*\/BulkRNASeq\/.*$/,
      },
      model_name: "rna_seq",
      watch_type: "link_files",
    )
  end
end
