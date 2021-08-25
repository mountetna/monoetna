require_relative "./ipi_rna_seq_files_linker_base"

class Polyphemus::IpiRnaSeqProcessedFilesLinker < Polyphemus::IpiRnaSeqFilesLinkerBase
  RECORD_NAME_REGEX = /.*\/(?<record_name>.*)\/.*$/

  CRAM_REGEX = /.*\.deduplicated\.cram$/
  CRAM_INDEX_REGEX = /.*\.deduplicated\.cram\.crai$/
  JUNCTION_REGEX = /.*\.junction$/
  UNMAPPED_FASTQ_REGEX = /.*\.fastq\.gz$/

  def initialize
    super(attribute_regex: attribute_regex, record_name_regex: RECORD_NAME_REGEX)
  end

  private

  def attribute_regex
    {
      "cram": CRAM_REGEX,
      "cram_index": CRAM_INDEX_REGEX,
      "junction": JUNCTION_REGEX,
      "unmapped_fastqs": UNMAPPED_FASTQ_REGEX,
    }
  end
end
