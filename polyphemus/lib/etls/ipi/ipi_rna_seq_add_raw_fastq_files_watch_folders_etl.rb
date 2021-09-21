require_relative "../add_watch_folder_base_etl"
require_relative "./linkers/ipi_rna_seq_raw_fastq_files_linker"

class Polyphemus::IpiRnaSeqAddRawFastqFilesWatchFoldersEtl < Polyphemus::AddWatchFolderBaseEtl
  PROJECT = "ipi"
  BUCKET = "integral_data"

  def initialize
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_path_regexes: {
        "#{PROJECT}_#{BUCKET}": /^.*\/BulkRNASeq\/.*$/,
      },
      model_name: "rna_seq",
      watch_type: "link_files",
    )
  end

  def should_skip_folder?(folder)
    # Ideally we wouldn't watch folders that aren't actual, valid record
    #   names in Magma. However, we are creating record names from the folders currently,
    #   so we would run into a race condition between the ETLs. For now, we'll treat all folders
    #   as watch folders and filter out records during the linking process.
    # !current_magma_record_names.include?(folder.folder_name)
    false
  end
end
