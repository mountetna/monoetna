require_relative "../../metis_file_in_watch_folder_etl"
require_relative "./file_processors/ipi_rna_seq_matrix_processor"

class Polyphemus::IpiRnaSeqPopulateMatricesEtl < Polyphemus::MetisFileInWatchFolderEtl
  PROJECT = "ipi"
  BUCKET = "data"

  def initialize
    @processor = Polyphemus::IpiRnaSeqMatrixProcessor.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      limit: 1,
      watch_type: "process_files",
      hide_paths: false,
    )
  end

  def process(cursor, files)
    logger.info("Found #{files.length} files in watch folders: #{files.map { |f| f.file_path }.join(",")}")

    @processor.process(cursor, files)

    logger.info("Done")
  end
end
