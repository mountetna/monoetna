require_relative "../link_files_base_etl"
require_relative "../../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqLinkProcessedFilesBaseEtl < Polyphemus::LinkFilesBaseEtl
  PLATE_REGEX = /.*\/(?<plate>plate\d+)_.*\/output\/.*/

  PROJECT = "ipi"
  BUCKET = "data"
  MODEL = "rna_seq"

  def initialize(attribute_name:, path_regex:)
    @helper = IpiHelper.new
    super(
      attribute_name: attribute_name,
      path_regex: path_regex,
      file_name_globs: ["bulkRNASeq/**/*", "output/**/*"],
      project_bucket_model_tuples: [[PROJECT, BUCKET, MODEL]],
    )
  end

  private

  def link_model_record_name(metis_file)
    metis_file.file_path.match(PLATE_REGEX)[:plate].gsub(/plate/, "Plate")
  end

  def corrected_record_name(link_model_record_name, record_name)
    @helper.corrected_rna_seq_tube_name(link_model_record_name, record_name)
  end

  def should_skip_record?(record_name)
    @helper.is_non_cancer_sample?(record_name)
  end
end
