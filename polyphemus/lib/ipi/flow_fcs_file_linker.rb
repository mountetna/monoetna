require 'etna'
require 'json'
require_relative './ipi_helper'

class IpiFlowFcsLinker < Etna::Clients::Magma::FileLinkingWorkflow
  def initialize(**opts)
    super(**{bucket_name: 'integral_data', model_name: 'flow'}.update(opts))
    @helper = IpiHelper.new
  end

  def find_matches
    {}.tap do |all_matches|
      @helper.find_files_by_extension(
        metis_client: metis_client,
        project_name: project_name,
        bucket_name: bucket_name,
        extension: 'fcs',
        flow_fcs_only: true
      ).each do |fcs_file|
        file_path = fcs_file.file_path
        key = [{
          fcs_file: file_path,
          stain: @helper.stain_from_fcs(file_path),
          stain_name: @helper.stain_name_from_fcs(file_path)
        }, :fcs_file]
        all_matches[key] = [file_path]
      end
    end
  end

  # Subclasses should override this to implement custom logic for how regex matches should match to linking.
  def matches_to_record_identifiers(match_data)
    super(match_data).tap do |identifiers|
      identifiers.update(match_data)
      # Below data is to get_or_create the parent records
      identifiers['flow'] = match_data[:stain_name]
      identifiers['sample'] = @helper.sample_name_from_fcs(match_data[:fcs_file])
      identifiers['patient'] = @helper.patient_name_from_fcs(match_data[:fcs_file])
      identifiers['experiment'] = @helper.experiment_from_patient_number(identifiers['patient'])
    end
  end

  def revision_for(id, attribute_name, file_path, match_map, record_identifiers)
    # Make sure the new flow record is linked correctly in the graph, and
    #   populate the stain attribute.
    super.update({
      sample: record_identifiers['sample'],
      stain: match_map[:stain]
    })
  end

  def link_files
    each_revision do |id, revision|
      # Not using magma_crud here with bulk updates, because was getting a
      #   vague Server Error when submitting all at once.
      # It might have been dev-environment / proxy related, but doing individual
      #   requests per revision for sanity's sake and to catch if it is an
      #   issue with Metis.
      # Have to inject project_name into the update request.
      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: 'ipi')
      update_request.update_revision(model_name, id, revision)
      puts "Updating #{model_name}, record #{id}."
      magma_crud.magma_client.update_json(update_request)
    end
  end
end