class IpiHelper
  def initialize
  end

  def find_files_by_extension(metis_client:, project_name:, bucket_name:, extension:)
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

  def find_files_in_folder_by_extension(metis_client:, project_name:, bucket_name:, extension:, folder_name:)
    metis_client.find(
      Etna::Clients::Metis::FindRequest.new(
        project_name: project_name,
        bucket_name: bucket_name,
        params: [Etna::Clients::Metis::FindParam.new(
          attribute: 'name',
          predicate: 'glob',
          value: "#{folder_name}/**/*.#{extension}",
          type: 'file'
        )]
    )).files.all
  end

  def experiment_from_patient_number(patient_number)
    experiments = JSON.parse(File.read('./lib/ipi/ipi_experiment_map.json'))
    exp_codes = experiments.keys
    reg = Regexp.new("(?<exp>#{exp_codes.join('|')})")
    experiments[reg.match(patient_number)[:exp]]
  end

  def file_name_from_path(path)
    path.split('/').last
  end

  def patient_name_from_wsp(wsp_file_path)
    file_name_from_path(wsp_file_path).sub('.wsp', '')
  end

  def stain_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return dc
    file_name = file_name_from_path(fcs_file_path).sub('.fcs', '')
    regex_results = fcs_regex(file_name)
    regex_results[:stain].downcase
  end

  def stain_name_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return IPICRC024.T3.flow.dc
    file_name_from_path(
      fcs_file_path).sub('.fcs', '').gsub('_', '.')
  end

  def sample_name_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return IPICRC024.T3
    file_name = file_name_from_path(fcs_file_path).sub('.fcs', '')
    regex_results = fcs_regex(file_name)
    regex_results[:sample_name].sub('_', '.')
  end

  def patient_name_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return IPICRC024
    sample_name = sample_name_from_fcs(fcs_file_path)
    sample_name.split('.').first
  end

  def fcs_regex(file_name)
    regex = Regexp.new('(?<sample_name>IPI\w+\d+_\w+\d)_(?:flow_)?(?<stain>dc|innate|nktb|sort|treg)', Regexp::IGNORECASE)
    raise "FCS filename #{file_name} does not match expected regex" unless regex.match(file_name)
    regex.match(file_name)
  end
end