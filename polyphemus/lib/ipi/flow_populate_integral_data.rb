# This class copies data from a Metis bucket and organizes it into the integral dataset.
# This assumes the WSP and FCS files in the source Metis bucket
#   have gone through Lenny's pipeline and are consistently named,
#   following our naming convention.
require 'ostruct'
require_relative './ipi_helper'

class IpiFlowPopulateIntegralData < Struct.new(:metis_client, :project_name, :source_bucket_name, :source_folder_name, :dest_bucket_name, keyword_init: true)
  attr_reader :helper

  def initialize(**opts)
    super(**{dest_bucket_name: 'integral_data', project_name: 'ipi'}.update(opts))
    @helper = IpiHelper.new
  end

  def find_files(extension)
    helper.find_files_in_folder_by_extension(
      metis_client: metis_client,
      project_name: project_name,
      bucket_name: source_bucket_name,
      folder_name: source_folder_name,
      extension: extension
    )
  end

  def source_metis_path(file)
    "metis://#{project_name}/#{source_bucket_name}/#{file.file_path}"
  end

  def dest_path(file)
    # Return the full metis path for the file
    #
    # WSP are associated with the patient
    #
    # IPIGYN001/
    #     IPIGYN001.wsp
    #
    # FCS are associated with the sample / tube, i.e. FlowCytometry/ would be
    #     a new directory under IPIGYN001.T1/
    # Then each stain has a subdirectory?
    #
    # IPIGYN001.T1/
    #     FlowCytometry/
    #         IPIGYN001.T1.flow.dc/
    #             IPIGYN001_T1_flow_dc.fcs

    patient_name = is_wsp(file) ?
      helper.patient_name_from_wsp(file.file_path) :
      helper.patient_name_from_fcs(file.file_path)
    experiment = helper.experiment_from_patient_number(patient_name)

    return "#{experiment}/#{patient_name}/#{file.file_name}" if is_wsp(file)

    sample_name = helper.sample_name_from_fcs(file.file_path)
    stain_name = helper.stain_name_from_fcs(file.file_path)

    return "#{experiment}/#{patient_name}/#{sample_name}/FlowCytometry/#{stain_name}/#{file.file_name}"
  end

  def dest_metis_path(file)
    return "metis://#{project_name}/#{dest_bucket_name}/#{dest_path(file)}"
  end

  def is_wsp(file)
      return file.file_name.downcase.end_with?('wsp')
  end

  def copy_files_to_dest(ext)
    find_files(ext).each do |file|
      copy_files_request = Etna::Clients::Metis::CopyFilesRequest.new(project_name: project_name)

      puts "Ensuring parent folder exists for #{file.file_name}."
      metis_client.ensure_parent_folder_exists(
        project_name: project_name,
        bucket_name: dest_bucket_name,
        path: dest_path(file))

      puts "Creating copy request."
      copy_files_request.add_revision(
        Etna::Clients::Metis::CopyRevision.new(
          source: source_metis_path(file),
          dest: dest_metis_path(file)
        )
      )

      metis_client.signed_copy_files(copy_files_request, Polyphemus.instance, 'metis')
    end
  end

  def copy_files
    ['wsp', 'fcs'].each do |ext|
      copy_files_to_dest(ext)
    end
  end
end