require_relative "./ipi_rna_seq_files_linker_base"

class Polyphemus::IpiRnaSeqRawFastqFilesLinker < Polyphemus::IpiRnaSeqFilesLinkerBase
  RECORD_NAME_REGEX = /.*\/(?<record_name>.*)\/.*$/

  RAW_FASTQ_FILE = /.*\.fastq\.gz$/

  def initialize
    super(attribute_regex: attribute_regex, record_name_regex: RECORD_NAME_REGEX)
  end

  private

  def attribute_regex
    {
      "raw_fastq_files": RAW_FASTQ_FILE,
    }
  end
end
