require_relative "./ipi_rna_seq_link_processed_files_base_etl"

class Polyphemus::IpiRnaSeqLinkUnmappedFastqFilesEtl < Polyphemus::IpiRnaSeqLinkProcessedFilesBaseEtl
  PATH_REGEX = /.*\/(?<record_name>.*)\/(?<original_file_name>.*\.fastq\.gz)$/

  def initialize
    super(
      attribute_name: "unmapped_fastqs",
      path_regex: PATH_REGEX,
    )
  end
end
