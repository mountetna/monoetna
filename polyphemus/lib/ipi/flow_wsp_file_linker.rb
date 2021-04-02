require 'etna'
require 'json'
require_relative './ipi_helper'

class IpiFlowWspLinker < Etna::Clients::Magma::FileLinkingWorkflow
  def initialize(**opts)
    super(**{bucket_name: 'integral_data', model_name: 'patient'}.update(opts))
    @helper = IpiHelper.new
  end

  def find_matches
    {}.tap do |all_matches|
      @helper.find_files_by_extension(
        metis_client: metis_client,
        project_name: project_name,
        bucket_name: bucket_name,
        extension: 'wsp'
      ).each do |wsp_file|
        file_path = wsp_file.file_path
        if file_path.include?("_ORIGINAL")
          key = [{
            ipi_number: @helper.patient_name_from_wsp(file_path),
            flojo_file_raw: file_path
          }, :flojo_file_raw]
        else
          key = [{
            ipi_number: @helper.patient_name_from_wsp(file_path),
            flojo_file_processed: file_path
          }, :flojo_file_processed]
        end
        all_matches[key] = [file_path]
      end
    end
  end

  # Subclasses should override this to implement custom logic for how regex matches should match to linking.
  def matches_to_record_identifiers(match_data)
    super(match_data).tap do |identifiers|
      identifiers.update(match_data)
      # Below data is to get_or_create the parent records
      identifiers['patient'] = match_data[:ipi_number]
      identifiers['experiment'] = @helper.experiment_from_patient_number(match_data[:ipi_number])
    end
  end

  def link_files
    each_revision do |id, revision|
      # Not using magma_crud here with bulk updates, because was getting a
      #   Net::ReadTimeout from Metis on the bulk_copy request.
      # Have to inject project_name into the update request.
      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: 'ipi')
      update_request.update_revision(model_name, id, revision)
      puts "Updating #{model_name}, record #{id}."
      magma_crud.magma_client.update(update_request)
    end
  end
end