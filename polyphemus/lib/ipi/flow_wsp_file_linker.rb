require 'etna'
require 'json'

class IpiFlowWspLinker < Etna::Clients::Magma::FileLinkingWorkflow
  def initialize(**opts)
    super(**{bucket_name: 'integral_data', model_name: 'patient'}.update(opts))
  end

  def find_matches
    {}.tap do |all_matches|
      metis_client.find(
        Etna::Clients::Metis::FindRequest.new(
          project_name: project_name,
          bucket_name: bucket_name,
          params: [Etna::Clients::Metis::FindParam.new(
            attribute: 'name',
            predicate: '=~',
            value: '/.*\.wsp$/i',
            type: 'file'
          )]
      )).files.all.each do |wsp_file|
        key = [{
          ipi_number: patient_name_from_wsp(wsp_file.file_path),
          flojo_file: wsp_file.file_path}, :flojo_file]
        all_matches[key] = [wsp_file.file_path]
      end
    end
  end

  def patient_name_from_wsp(wsp_file_path)
    wsp_file_path.split('/').last.sub('.wsp', '')
  end

  def experiment_from_patient_number(patient_number)
    experiments = JSON.parse(File.read('./lib/ipi/ipi_experiment_map.json'))
    exp_codes = experiments.keys
    reg = Regexp.new("(?<exp>#{exp_codes.join('|')})")
    experiments[reg.match(patient_number)[:exp]]
  end

  # Subclasses should override this to implement custom logic for how regex matches should match to linking.
  def matches_to_record_identifiers(match_data)
    super(match_data).tap do |identifiers|
      identifiers.update(match_data)
      # Below data is to get_or_create the parent records
      identifiers['patient'] = match_data[:ipi_number]
      identifiers['experiment'] = experiment_from_patient_number(match_data[:ipi_number])
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