require_relative "../add_watch_folder_base_etl"
require_relative "./ipi_rna_seq_processed_files_linker"
require_relative "../../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqPopulateProcessedFilesWatchFoldersEtl < Polyphemus::AddWatchFolderBaseEtl
  PROJECT = "ipi"
  BUCKET = "data"
  MODEL = "rna_seq"

  def initialize
    @linker = Polyphemus::IpiRnaSeqProcessedFilesLinker.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_name_globs: ["output/*", "bulkRNASeq/*"],
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

  def link_folder_contents(cursor, folders)
    folders.each do |folder|
      files = list_folder_files(cursor, folder)
      @linker.link(
        project_name: cursor[:project_name],
        model_name: MODEL,
        files: files,
      )
    end
  end

  def list_folder_files(cursor, folder)
    self.metis_client.list_folder(Etna::Clients::Metis::ListFolderRequest.new(
      project_name: cursor[:project_name],
      bucket_name: cursor[:bucket_name],
      folder_path: folder.folder_path,
    )).files.all
  end
end
