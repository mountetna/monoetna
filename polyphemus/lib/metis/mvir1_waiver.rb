require 'etna'
require_relative '../commands'

class Mvir1Waiver
  include Polyphemus::WithEtnaClients
  include Polyphemus::WithLogger

  def initialize(metis_client:, project_name: 'mvir1')
    @metis_client = metis_client
    @project_name = project_name
  end

  def patient_name_regex(patient_name)
    # This expects a Patient name like MVIR1-HS10,
    #   not an assay name, like MVIR1-HS10-D4BLD1-CYM2
    /.*\/#{patient_name}-[\w-]+$/
  end

  def pool_regex(pool_name)
    # Slightly different regex for pools because we should be able
    #   to get the whole pool_name instead of a patient_name
    #   that we can't match to folder names.
    /.*\/#{pool_name}$/
  end

  def restrict_patient_data(patient_name)
    rename_folders_by_regex(
      source_bucket: bucket_name(:release_bucket),
      source_folders: folders(:release_bucket),
      dest_bucket: bucket_name(:restrict_bucket),
      regex: patient_name_regex(patient_name))
  end

  def release_patient_data(patient_name)
    rename_folders_by_regex(
      source_bucket: bucket_name(:restrict_bucket),
      source_folders: folders(:restrict_bucket),
      dest_bucket: bucket_name(:release_bucket),
      regex: patient_name_regex(patient_name))
  end

  def restrict_pool_data(pool_name)
    rename_folders_by_regex(
      source_bucket: bucket_name(:release_bucket),
      source_folders: folders(:release_bucket),
      dest_bucket: bucket_name(:restrict_bucket),
      regex: pool_regex(pool_name))
  end

  def release_pool_data(pool_name)
    rename_folders_by_regex(
      source_bucket: bucket_name(:restrict_bucket),
      source_folders: folders(:restrict_bucket),
      dest_bucket: bucket_name(:release_bucket),
      regex: pool_regex(pool_name))
  end

  private

  def bucket_name(bucket_name)
    Polyphemus.instance.config(:metis)[bucket_name]
  end

  def folders(bucket_name)
    @folders ||= Hash.new { |h, key|
      h[key] = @metis_client.list_all_folders(
        Etna::Clients::Metis::ListFoldersRequest.new(
          project_name: @project_name,
          bucket_name: bucket_name(key)
        )).folders.all
    }

    @folders[bucket_name]
  end

  def rename_folders_by_regex(source_bucket:, source_folders:, dest_bucket:, regex:)
    found_folders = source_folders.select { |folder|
        folder.folder_path =~ regex
    }

    msg = "No found folders in #{source_bucket} matching #{regex} -- was a valid folder name provided?"
    return logger.warn(msg) if found_folders.length == 0

    found_folders.each { |folder|
        @metis_client.rename_folder(Etna::Clients::Metis::RenameFolderRequest.new(
          bucket_name: source_bucket,
          project_name: @project_name,
          folder_path: folder.folder_path,
          new_bucket_name: dest_bucket,
          new_folder_path: folder.folder_path,
          create_parent: true)
        )
    }
  end
end