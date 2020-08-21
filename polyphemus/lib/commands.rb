require 'date'
require 'logger'

require 'pry'

class Polyphemus
  class Help < Etna::Command
    usage 'List this help'

    def execute
      puts 'Commands:'
      Polyphemus.instance.commands.each do |name,cmd|
        puts cmd.usage
      end
    end
  end

  module WithEtnaClients
    def project
      raise "project must be implemented in subclasses!"
    end

    def token
      Polyphemus.instance.config(project)[:token]
    end

    def magma_client
      @magma_client ||= Etna::Clients::Magma.new(token: token, host: Polyphemus.instance.config(:magma)[:host])
    end

    def metis_client
      @metis_client ||= Etna::Clients::Metis.new(token: token, host: Polyphemus.instance.config(:metis)[:host])
    end
  end

  module WithLogger
    def logger
      Polyphemus.instance.logger
    end
  end

  module WithMetisWaiverHelpers
    include WithEtnaClients
    include WithLogger

    def project
      raise "project must be implemented in subclasses!"
    end

    def restrict_patient_data(patient_name)
      # Move a folder from release_bucket to restrict_bucket
      # If can't find the patient_name folder, throw exception.
      # If we use paging, this would also have to fetch more pages.
      # This expects a Patient name like MVIR1-HS10,
      #   not an assay name, like MVIR1-HS10-D4BLD1-CYM2
      patient_folders = release_folders.select { |folder|
          folder.folder_path =~ /.*\/#{patient_name}-\w+$/
      }

      # If we want to check by page, we would need to do
      #   fetch more pages here. However, have to
      #   make sure that we find the patient_name, so basically
      #   have to fetch all folders...
      msg = "No found folders in #{release_bucket_name} -- is #{patient_name} a valid patient name?"
      logger.error(msg) if patient_folders.length == 0
      raise Etna::Error, msg if patient_folders.length == 0

      patient_folders.each { |folder|
          rename_folder(
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
          folder.folder_path =~ /.*\/#{patient_name}-\w+$/
      }

      # If we want to check by page, we would need to do
      #   fetch more pages here. However, have to
      #   make sure that we find the patient_name, so basically
      #   have to fetch all folders...
      msg = "No found folders in #{restrict_bucket_name} -- is #{patient_name} a valid patient name?"
      logger.error(msg) if patient_folders.length == 0
      raise Etna::Error, msg if patient_folders.length == 0

      patient_folders.each { |folder|
          rename_folder(
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
          folder.folder_path =~ /.*\/#{pool_name}$/
      }

      # If we want to check by page, we would need to do
      #   fetch more pages here. However, have to
      #   make sure that we find the patient_name, so basically
      #   have to fetch all folders...
      msg = "No found folders in #{release_bucket_name} -- is #{pool_name} a valid pool name?"
      logger.error(msg) if pool_folders.length == 0
      raise Etna::Error, msg if pool_folders.length == 0

      pool_folders.each { |folder|
          rename_folder(
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
          folder.folder_path =~ /.*\/#{pool_name}$/
      }

      # If we want to check by page, we would need to do
      #   fetch more pages here. However, have to
      #   make sure that we find the patient_name, so basically
      #   have to fetch all folders...
      msg = "No found folders in #{restrict_bucket_name} -- is #{pool_name} a valid pool name?"
      logger.error(msg) if pool_folders.length == 0
      raise Etna::Error, msg if pool_folders.length == 0

      pool_folders.each { |folder|
          rename_folder(
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
      @release_folders ||= fetch_folders(release_bucket_name).all
    end

    def restrict_folders
      @restrict_folders ||= fetch_folders(restrict_bucket_name).all
    end

    def fetch_folders(bucket_name)
      metis_client.list_all_folders(
        Etna::Clients::Metis::ListFoldersRequest.new(
          project_name: project, bucket_name: bucket_name)).folders
    end

    def parent_folder_path(folder_path)
      folder_path.split('/')[0..-2].join('/')
    end

    def create_parent_folder(bucket_name, folder_path)
      parent_path = parent_folder_path(folder_path)
      logger.debug("Creating parent folder #{bucket_name}/#{folder_path}")
      metis_client.create_folder(
        Etna::Clients::Metis::CreateFolderRequest.new(
          project_name: project,
          bucket_name: bucket_name,
          folder_path: parent_path
      ))
    end

    def parent_exists?(bucket_name, folder_path)
      # NOTE: this doesn't test if the folder_path itself exists
      #   This can be confusing for root folders, because
      #       they have no parents, so you don't need
      #       to create anything.
      parent_path = parent_folder_path(folder_path)
      return true if parent_path.empty? # root folder

      # returns 422 if the folder_path does not exist
      begin
          metis_client.list_folder(
            Etna::Clients::Metis::ListFolderRequest.new(
              project_name: project,
              bucket_name: bucket_name,
              folder_path: parent_path
          ))
      rescue Etna::Error => e
          return false if e.status == 422
          raise
      end
      return true
    end

    def rename_folder(source_bucket_name, source_folder_path, dest_bucket_name, dest_folder_path)
      logger.debug("Renaming folder #{source_bucket_name}/#{source_folder_path} to #{dest_bucket_name}/#{dest_folder_path}")
      create_parent_folder(dest_bucket_name, dest_folder_path) if !parent_exists?(dest_bucket_name, dest_folder_path)

      logger.debug("Making the call to folder_rename route")
      metis_client.rename_folder(
        Etna::Clients::Metis::RenameFolderRequest.new(
          bucket_name: source_bucket_name,
          project_name: project,
          folder_path: source_folder_path,
          new_bucket_name: dest_bucket_name,
          new_folder_path: dest_folder_path
      ))
    end
  end

  class CascadeMvirPatientWaiverToRestricted < Etna::Command
    include WithEtnaClients
    include WithLogger
    include WithMetisWaiverHelpers

    usage 'Updates any models whose restricted access does not match its relationship to patient'

    def project
      :mvir1
    end

    def execute
      request = Etna::Clients::Magma::RetrievalRequest.new(project_name: project)
      request.model_name = 'patient'
      request.attribute_names = 'all'
      request.record_names = 'all'
      patient_documents = magma_client.retrieve(request).models.model('patient').documents

      # Each patient is handled so that the metis interactions can fire for individual patients and their status
      # updated upon successful individually.  This is technically less efficient than, say, batching one very large
      # update, but it has a crucial benefit in that any failures will leave the system in partial completion and thus
      # reduce the work a future run would need to complete as a result of individual patient states being aligned.
      patient_documents.document_keys.each do |patient_name|
        patient = patient_documents.document(patient_name)
        should_be_restricted = patient['consent'] == 'Initial Waiver'
        should_be_deleted = patient['consent'] == 'No Samples/Data'

        if should_be_deleted
          # Magma does not support delete right now, so we restrict the models, but we can unlink and
          # delete the files in metis
          restrict!(patient, delete_metis_files: true)
        elsif should_be_restricted && ! patient['restricted']
          restrict!(patient)
        elsif ! should_be_restricted && patient['restricted']
          unrestrict!(patient)
        end
      end

      cascade_to_pools
    end

    def cascade_to_pools
      request = Etna::Clients::Magma::QueryRequest.new(project_name: project)
      request.query = [ 'cytof',
                        [ 'timepoint', 'patient', 'restricted', '::true' ],
                        '::all', 'cytof_pool', '::identifier' ]
      all_restricted_pools = magma_client.query(request).answer.map { |r| r[1] }.sort.uniq

      request.query = [ 'cytof_pool', '::all', '::identifier' ]
      all_pools = magma_client.query(request).answer.map { |r| r[1] }.sort.uniq

      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)

      all_pools.each do |pool|
        if all_restricted_pools.include? pool
          logger.info "Cytof pool #{pool} includes a restricted patient, restricting."

          # This should be done per-pool, like patients, not bulk?
          restrict_pool_data(pool)

          update_request.update_revision('cytof_pool', pool, restricted: true)
        else
          logger.info "Cytof pool #{pool} does not include a restricted patient, relaxing."

          # This should be done per-pool, like patients, not bulk?
          release_pool_data(pool)

          update_request.update_revision('cytof_pool', pool, restricted: false)
        end
      end

      magma_client.update(update_request)
    end

    def restrict!(patient, delete_metis_files: false)
      name = patient['name']
      logger.warn("Attempting to restrict access to #{name}")

      # This code path should be --eventually consistent--  That is to say, we should ensure each operation
      # is idempotent (may need to be repeated), and that the patient is not marked restricted until all other
      # related tasks are complete and the state is consistent.
      # Metis update here!

      if delete_metis_files
        # do metis deletion attempt
        restrict_patient_data(name)
      else
        # do metis movement attempt
        restrict_patient_data(name)
      end


      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
      update_request.update_revision('patient', name, restricted: true)
      magma_client.update(update_request)
    end

    def unrestrict!(patient)
      name = patient['name']
      logger.warn("Attempting to unrestrict access to #{name}")

      # This code path should be --eventually consistent--  That is to say, we should ensure each operation
      # is idempotent (may need to be repeated), and that the patient is not marked restricted until all other
      # related tasks are complete and the state is consistent.
      # Metis update here!
      release_patient_data(name)

      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
      update_request.update_revision('patient', name, restricted: false)
      magma_client.update(update_request)
    end

    def setup(config)
      super
      Polyphemus.instance.setup_logger
    end
  end

  class GetMetisFolders < Etna::Command
    include WithEtnaClients

    usage "Fetch a list of Metis folders from a bucket"

    def project
      :mvir1
    end

    def execute
      folders = metis_client.list_all_folders(
        Etna::Clients::Metis::ListFoldersRequest.new(project_name: project, bucket_name: 'data')).folders
      folders.all.each { |f| p f.folder_path }
    end

    def setup(config)
      super
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Polyphemus instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
    end
  end
end
