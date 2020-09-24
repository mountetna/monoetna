require 'etna'
require 'json'
require_relative './ipi_helper'

class IpiFlowWspLinker < Etna::Clients::Magma::FileLinkingWorkflow
  def initialize(**opts)
    super(**{bucket_name: 'integral_data', model_name: 'patient'}.update(opts))
  end

  def find_matches
    {}.tap do |all_matches|
      IpiHelper::find_files_by_extension(
        metis_client: metis_client,
        project_name: project_name,
        bucket_name: bucket_name,
        extension: 'wsp'
      ).each do |wsp_file|
        file_path = wsp_file.file_path
        key = [{
          ipi_number: patient_name_from_wsp(file_path),
          flojo_file: file_path
        }, :flojo_file]
        all_matches[key] = [file_path]
      end
    end
  end

  def patient_name_from_wsp(wsp_file_path)
    IpiHelper::file_name_from_path(wsp_file_path).sub('.wsp', '')
  end

  # Subclasses should override this to implement custom logic for how regex matches should match to linking.
  def matches_to_record_identifiers(match_data)
    super(match_data).tap do |identifiers|
      identifiers.update(match_data)
      # Below data is to get_or_create the parent records
      identifiers['patient'] = match_data[:ipi_number]
      identifiers['experiment'] = IpiHelper::experiment_from_patient_number(match_data[:ipi_number])
    end
  end

  def link_files
    magma_crud.update_records do |update_request|
      each_revision do |id, revision|
        update_request.update_revision(model_name, id, revision)
      end
    end
  end
end