require_relative "../add_watch_folder_base_etl"
require_relative "./file_processors/ipi_rna_seq_attribute_processor"
require_relative "./file_processors/ipi_rna_seq_matrix_processor"

class Polyphemus::IpiRnaSeqAddResultsWatchFoldersEtl < Polyphemus::AddWatchFolderBaseEtl
  PROJECT = "ipi"
  BUCKET = "data"

  def initialize
    @attribute_processor = Polyphemus::IpiRnaSeqAttributeProcessor.new
    @matrix_processor = Polyphemus::IpiRnaSeqMatrixProcessor.new
    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_path_regexes: {
        "#{PROJECT}_#{BUCKET}": /^bulkRNASeq\/.*\/results$/,
      },
      model_name: "rna_seq",
      watch_type: "process_files",
    )
  end

  def process_folder_contents(cursor, folders)
    # Need to actually download and process several TSV files,
    #   so override and don't just link here.
    folders.each do |folder|
      files = folder_files(cursor, folder)
      @attribute_processor.process(cursor, files)
      @matrix_processor.process(cursor, files)
    end
  end
end
