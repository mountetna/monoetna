require_relative "../metis_files_linker_base"
require_relative "../ipi/ipi_helper"
require_relative "../../helpers"

class Polyphemus::IpiRnaSeqProcessedFilesLinker < Polyphemus::MetisFilesLinkerBase
  WithEtnaClients

  PATH_REGEX = /.*\/(?<record_name>.*)\/.*$/

  CRAM_REGEX = /.*\.deduplicated\.cram$/
  CRAM_INDEX_REGEX = /.*\.deduplicated\.cram\.crai$/
  JUNCTION_REGEX = /.*\.junction$/
  UNMAPPED_FASTQ_REGEX = /.*\.fastq\.gz$/

  def initialize
    super
    @helper = IpiHelper.new
  end

  def link(project_name:, model_name:, files:)
    super(
      project_name: project_name,
      model_name: model_name,
      files_by_record_name: Polyphemus::MetisFilesLinkerBase.organize_metis_files_by_magma_record(
        metis_files: files,
        magma_record_names: current_magma_record_names(project_name, model_name),
        path_regex: PATH_REGEX,
      ),
      attribute_regex: attribute_regex,
    )
  end

  def corrected_record_name(record_name)
    @helper.corrected_rna_seq_tube_name(record_name)
  end

  def should_skip_record?(record_name)
    @helper.is_non_cancer_sample?(record_name)
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

  def current_magma_record_names(project_name, model_name)
    @current_magma_record_names ||= begin
        request = Etna::Clients::Magma::RetrievalRequest.new(
          project_name: project_name,
          model_name: model_name,
          attribute_names: ["identifier"],
          record_names: "all",
          hide_templates: true,
        )
        self.magma_client.retrieve(request).models.model(model_name).documents.document_keys
      end
  end
end
