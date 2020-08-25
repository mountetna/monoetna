require 'etna'
require_relative '../commands'

class MetisWaiverHelper
  include Polyphemus::WithEtnaClients
  include Polyphemus::WithLogger

  def initialize(metis_client)
    @metis_client = metis_client
  end

  def patient_name_regex(patient_name)
    /.*\/#{patient_name}-[\w-]+$/
  end

  def pool_regex(pool_name)
    /.*\/#{pool_name}$/
  end

  def restrict_patient_data(patient_name)
    # Move a folder from release_bucket to restrict_bucket
    # If can't find the patient_name folder, throw exception.
    # If we use paging, this would also have to fetch more pages.
    # This expects a Patient name like MVIR1-HS10,
    #   not an assay name, like MVIR1-HS10-D4BLD1-CYM2
    patient_folders = release_folders.select { |folder|
        folder.folder_path =~ patient_name_regex(patient_name)
    }

    # If we want to check by page, we would need to do
    #   fetch more pages here. However, have to
    #   make sure that we find the patient_name, so basically
    #   have to fetch all folders...
    msg = "No found folders in #{release_bucket_name} -- is #{patient_name} a valid patient name?"
    return logger.warn(msg) if patient_folders.length == 0

    patient_folders.each { |folder|
        @metis_client.rename_folder(
            release_bucket_name,
            folder.folder_path,
            restrict_bucket_name,
            folder.folder_path
        )
    }
  end

  def release_patient_data(patient_name)
    # Move a folder from restrict_bucket to release_bucket
    # If can't find the patient_name folder, throw exception.
    # If we use paging, this would also have to fetch more pages.
    # This expects a Patient name like MVIR1-HS10,
    #   not an assay name, like MVIR1-HS10-D4BLD1-CYM2
    patient_folders = restrict_folders.select { |folder|
        folder.folder_path =~ patient_name_regex(patient_name)
    }

    # If we want to check by page, we would need to do
    #   fetch more pages here. However, have to
    #   make sure that we find the patient_name, so basically
    #   have to fetch all folders...
    msg = "No found folders in #{restrict_bucket_name} -- is #{patient_name} a valid patient name?"
    return logger.warn(msg) if patient_folders.length == 0

    patient_folders.each { |folder|
        @metis_client.rename_folder(
            restrict_bucket_name,
            folder.folder_path,
            release_bucket_name,
            folder.folder_path
        )
    }
  end

  def restrict_pool_data(pool_name)
    # Move a folder from release_bucket to restrict_bucket
    # Slightly different regex for pools because we should be able
    #   to get the whole pool_name instead of a patient_name
    #   that we can't match to folder names.
    pool_folders = release_folders.select { |folder|
        folder.folder_path =~ pool_regex(pool_name)
    }

    # If we want to check by page, we would need to do
    #   fetch more pages here. However, have to
    #   make sure that we find the patient_name, so basically
    #   have to fetch all folders...
    msg = "No found folders in #{release_bucket_name} -- is #{pool_name} a valid pool name?"
    return logger.warn(msg) if pool_folders.length == 0

    pool_folders.each { |folder|
        @metis_client.rename_folder(
            release_bucket_name,
            folder.folder_path,
            restrict_bucket_name,
            folder.folder_path
        )
    }
  end

  def release_pool_data(pool_name)
    # Move a folder from restrict_bucket to release_bucket
    # Slightly different regex for pools because we should be able
    #   to get the whole pool_name instead of a patient_name
    #   that we can't match to folder names.
    pool_folders = restrict_folders.select { |folder|
        folder.folder_path =~ pool_regex(pool_name)
    }

    # If we want to check by page, we would need to do
    #   fetch more pages here. However, have to
    #   make sure that we find the patient_name, so basically
    #   have to fetch all folders...
    msg = "No found folders in #{restrict_bucket_name} -- is #{pool_name} a valid pool name?"
    return logger.warn(msg) if pool_folders.length == 0

    pool_folders.each { |folder|
        @metis_client.rename_folder(
            restrict_bucket_name,
            folder.folder_path,
            release_bucket_name,
            folder.folder_path
        )
    }
  end

  private

  def release_bucket_name
    @release_bucket_name ||= Polyphemus.instance.config(:metis)[:release_bucket]
  end

  def restrict_bucket_name
    @restrict_bucket_name ||= Polyphemus.instance.config(:metis)[:restrict_bucket]
  end

  def release_folders
    @release_folders ||= @metis_client.fetch_folders(release_bucket_name).all
  end

  def restrict_folders
    @restrict_folders ||= @metis_client.fetch_folders(restrict_bucket_name).all
  end
end