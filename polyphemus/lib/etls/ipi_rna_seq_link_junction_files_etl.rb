require_relative "./ipi_rna_seq_link_processed_files_base_etl"

class Polyphemus::IpiRnaSeqLinkJunctionFilesEtl < Polyphemus::IpiRnaSeqLinkProcessedFilesBaseEtl
  PATH_REGEX = /.*\/(?<record_name>.*)\/(?<original_file_name>.*\.junction)$/

  def initialize
    super(
      attribute_name: "junction",
      path_regex: PATH_REGEX,
    )
  end
end
