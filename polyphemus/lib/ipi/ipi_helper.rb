require_relative "../data_processing/file_name_regex_dsl"

class IpiHelper
  include FileNameRegexDsl

  def initialize(renaming_file = "./lib/etls/renaming/projects/ipi_bulk_rna_renames.json")
    @renaming_file = renaming_file
  end

  def find_files_by_extension(metis_client:, project_name:, bucket_name:, extension:, flow_fcs_only: false)
    regex = "/.*\.#{extension}$/i" unless flow_fcs_only
    regex = "/.*_flow_.*\.#{extension}$/i" if flow_fcs_only

    metis_client.find(
      Etna::Clients::Metis::FindRequest.new(
        project_name: project_name,
        bucket_name: bucket_name,
        params: [Etna::Clients::Metis::FindParam.new(
          attribute: "name",
          predicate: "=~",
          value: regex,
          type: "file",
        )],
      )
    ).files.all
  end

  def find_files_in_folder_by_extension(metis_client:, project_name:, bucket_name:, extension:, folder_name:)
    metis_client.find(
      Etna::Clients::Metis::FindRequest.new(
        project_name: project_name,
        bucket_name: bucket_name,
        params: [Etna::Clients::Metis::FindParam.new(
          attribute: "name",
          predicate: "glob",
          value: "#{folder_name}/**/*.#{extension}",
          type: "file",
        )],
      )
    ).files.all
  end

  def experiment_from_patient_number(patient_number)
    experiments = JSON.parse(File.read("./lib/ipi/ipi_experiment_map.json"))
    exp_codes = experiments.keys
    reg = Regexp.new("(?<exp>#{exp_codes.join("|")})")
    experiments[reg.match(patient_number)[:exp]]
  end

  def file_name_from_path(path)
    path.split("/").last
  end

  def patient_name_from_wsp(wsp_file_path)
    file_name_from_path(wsp_file_path).sub(".wsp", "").sub("_PROCESSED", "")
  end

  def stain_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return dc
    file_name = file_name_from_path(fcs_file_path).sub(".fcs", "")
    regex_results = fcs_regex(file_name)
    regex_results[:stain].downcase
  end

  def stain_name_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return IPICRC024.T3.flow.dc
    file_name_from_path(
      fcs_file_path
    ).sub(".fcs", "").gsub("_", ".")
  end

  def sample_name_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return IPICRC024.T3
    file_name = file_name_from_path(fcs_file_path).sub(".fcs", "")
    regex_results = fcs_regex(file_name)
    regex_results[:sample_name].sub("_", ".")
  end

  def patient_name_from_fcs(fcs_file_path)
    # File name should be: IPICRC024_T3_flow_dc.fcs
    # We need to return IPICRC024
    sample_name = sample_name_from_fcs(fcs_file_path)
    sample_name.split(".").first
  end

  STAIN_NAMES = ["dc", "innate", "nktb", "sort", "treg"]

  def fcs_regex(file_name)
    regex = regex_chain(
      regex_named_match("sample_name", sample_name_regex("IPI")),
      /(?:flow_)?#{regex_named_match("stain", regex_options(STAIN_NAMES)).source}/,
      Regexp::IGNORECASE
    )
    raise "FCS filename #{file_name} does not match expected regex" unless regex.match(file_name)
    regex.match(file_name)
  end

  def corrected_rna_seq_tube_name(plate_name, tube_name)
    return tube_name unless rna_seq_renames.key?(plate_name)

    plate_renames = rna_seq_renames[plate_name]

    return tube_name unless plate_renames.key?(tube_name)

    plate_renames[tube_name]
  end

  def incorrect_rna_seq_tube_name(plate_name, tube_name)
    return tube_name unless rna_seq_renames.key?(plate_name)

    plate_renames_inverted = rna_seq_renames[plate_name].invert

    return tube_name unless plate_renames_inverted.key?(tube_name)

    plate_renames_inverted[tube_name]
  end

  def rna_seq_renames
    return {} unless ::File.exist?(@renaming_file)
    
    @rna_seq_renames ||= JSON.parse(File.read(@renaming_file))
  end

  def is_control?(record_name)
    record_name =~ /^control/i
  end

  def control_name(folder_name)
    # Control_(UHR|Jurkat).Plate\d+
    control, plate = folder_name.split(".")
    _, control_type = control.split("_")

    "Control_#{control_type =~ /jurkat/i ? "Jurkat" : "UHR"}.#{plate.capitalize}"
  end
end
