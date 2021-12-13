require_relative "../../../metis_files_linker_base"

class Polyphemus::RnaSeqFilesLinkerBase < Polyphemus::MetisFilesLinkerBase
  def initialize(project_name:, bucket_name:, attribute_regex:, record_name_regex:)
    super(project_name: project_name, bucket_name: bucket_name, record_name_regex: record_name_regex)
    @attribute_regex = attribute_regex
  end

  def link(model_name:, files:, cursor: nil)
    super(
      model_name: model_name,
      files_by_record_name: organize_metis_files_by_magma_record(
        metis_files: files,
        magma_record_names: current_magma_record_names(project_name, model_name),
      ),
      attribute_regex: @attribute_regex,
      cursor: cursor,
    )
  end
end
