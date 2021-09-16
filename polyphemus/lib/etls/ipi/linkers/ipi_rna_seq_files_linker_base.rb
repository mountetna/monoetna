require_relative "../../metis_files_linker_base"
require_relative "../../../ipi/ipi_helper"

class Polyphemus::IpiRnaSeqFilesLinkerBase < Polyphemus::MetisFilesLinkerBase
  def initialize(project_name:, bucket_name:, attribute_regex:, record_name_regex:)
    super(project_name: project_name, bucket_name: bucket_name)
    @attribute_regex = attribute_regex
    @record_name_regex = record_name_regex
    @helper = IpiHelper.new
  end

  def link(model_name:, files:)
    super(
      model_name: model_name,
      files_by_record_name: organize_metis_files_by_magma_record(
        model_name: model_name,
        metis_files: files,
        magma_record_names: current_magma_record_names(project_name, model_name),
        path_regex: @record_name_regex,
        attribute_regex: @attribute_regex,
      ),
      attribute_regex: @attribute_regex,
    )
  end

  private

  def corrected_record_name(record_name)
    @helper.corrected_rna_seq_tube_name(record_name)
  end

  def should_skip_record?(record_name)
    @helper.is_non_cancer_sample?(record_name)
  end

  def current_magma_record_names(project_name, model_name)
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
