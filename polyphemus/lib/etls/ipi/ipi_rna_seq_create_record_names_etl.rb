require_relative "../metis_folder_filtering_base_etl"
require_relative "./record_creators/ipi_rna_seq_and_plate_record_creator"

class Polyphemus::IpiRnaSeqCreateRecordNamesEtl < Polyphemus::MetisFolderFilteringBaseEtl
  PROJECT = "ipi"
  BUCKET = "data"

  def initialize
    @record_creator = Polyphemus::IpiRnaSeqAndPlateRecordCreator.new

    super(
      project_bucket_pairs: [[PROJECT, BUCKET]],
      folder_path_regexes: {
        "#{PROJECT}_#{BUCKET}": /^bulkRNASeq\/.*\/output\/.*$/,
      },
    )
  end

  def process(cursor, folders)
    # We'll need to filter out the folders based on folder_name here, using the
    #   supplied regexes
    logger.info("Found #{folders.length} updated folders: #{folders.map { |f| f.folder_path }.join(",")}")

    target_folders = filter_target_folders(cursor, folders)

    logger.info("Found #{target_folders.length} folders that match the cursor regex: #{target_folders.map { |f| f.folder_path }.join(",")}")

    process_folders(target_folders) unless target_folders.empty?

    logger.info("Done")
  end

  private

  def process_folders(folders)
    @record_creator.create(folders)
  end
end
