require_relative "../../metis_file_in_watch_folder_etl"
require_relative "./linkers/ipi_rna_seq_raw_fastq_files_linker"

class Polyphemus::IpiRnaSeqLinkRawFastqFilesEtl < Polyphemus::MetisFileInWatchFolderEtl
  PROJECT = "ipi"
  BUCKET = "integral_data"
  MODEL = "rna_seq"

  def initialize
    @linker = Polyphemus::IpiRnaSeqRawFastqFilesLinker.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      limit: 200,
      file_name_params: {
        "=~": ["*.fastq.gz"]
      }
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
