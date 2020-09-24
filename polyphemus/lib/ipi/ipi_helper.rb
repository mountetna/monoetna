class IpiHelper
  def self.find_files_by_extension(metis_client:, project_name:, bucket_name:, extension:)
    metis_client.find(
      Etna::Clients::Metis::FindRequest.new(
        project_name: project_name,
        bucket_name: bucket_name,
        params: [Etna::Clients::Metis::FindParam.new(
          attribute: 'name',
          predicate: '=~',
          value: "/.*\.#{extension}$/i",
          type: 'file'
        )]
    )).files.all
  end

  def self.experiment_from_patient_number(patient_number)
    experiments = JSON.parse(File.read('./lib/ipi/ipi_experiment_map.json'))
    exp_codes = experiments.keys
    reg = Regexp.new("(?<exp>#{exp_codes.join('|')})")
    experiments[reg.match(patient_number)[:exp]]
  end

  def self.file_name_from_path(path)
    path.split('/').last
  end
end