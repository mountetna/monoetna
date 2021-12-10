require_relative "../../shared/bulk_rna_seq/linkers/rna_seq_files_linker_base"
require_relative "../../../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqFilesLinkerBase < Polyphemus::RnaSeqFilesLinkerBase
  def initialize(project_name: "ipi", bucket_name:, attribute_regex:, record_name_regex:)
    super
    @helper = IpiHelper.new
  end

  def corrected_record_name(record_name)
    @helper.corrected_rna_seq_tube_name(record_name)
  end

  def should_skip_record?(record_name)
    @helper.is_non_cancer_sample?(record_name)
  end
end
