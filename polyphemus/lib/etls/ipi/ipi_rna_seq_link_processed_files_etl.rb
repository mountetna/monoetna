require_relative "../../metis_file_in_watch_folder_etl"
require_relative "./linkers/ipi_rna_seq_processed_files_linker"

class Polyphemus::IpiRnaSeqLinkProcessedFilesEtl < Polyphemus::MetisFileInWatchFolderEtl
  PROJECT = "ipi"
  BUCKET = "data"
  MODEL = "rna_seq"

  def initialize
    @linker = Polyphemus::IpiRnaSeqProcessedFilesLinker.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      watch_type: "link_files",
      hide_paths: true,
    )
  end

  def process(cursor, files)
    logger.info("Found #{files.length} files in watch folders: #{files.map { |f| f.file_path }.join(",")}")

    @linker.link(
      project_name: cursor[:project_name],
      model_name: MODEL,
      files: files,
    )

    logger.info("Done")
  end
end
