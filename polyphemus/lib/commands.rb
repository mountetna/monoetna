require 'date'
require 'logger'
require 'rollbar'


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
      Polyphemus.instance.config(:polyphemus)[:token]
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

  class CascadeMvirPatientWaiverToRestricted < Etna::Command
    include WithEtnaClients
    include WithLogger

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
        elsif should_be_restricted
          restrict!(patient)
        elsif ! should_be_restricted && patient['restricted']
          unrestrict!(patient)
        end
      end

      # Check all pools and verify status
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
      all_pools.each do |pool|
        if all_restricted_pools.include? pool
          logger.info "Cytof pool #{pool} includes a restricted patient, restricting."

          mvir1_waiver.restrict_pool_data(pool)

          update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
          update_request.update_revision('cytof_pool', pool, restricted: true)
          magma_client.update(update_request)
        else
          logger.info "Cytof pool #{pool} does not include a restricted patient, relaxing."

          mvir1_waiver.release_pool_data(pool)

          update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
          update_request.update_revision('cytof_pool', pool, restricted: false)
          magma_client.update(update_request)
        end
      end
    end

    def restrict!(patient, delete_metis_files: false)
      name = patient['name']

      if delete_metis_files
        # do metis movement attempt for now -- no auto-delete
        mvir1_waiver.restrict_patient_data(name)
      else
        # do metis movement attempt
        mvir1_waiver.restrict_patient_data(name)
      end

      unless patient['restricted']
        logger.warn("Attempting to restrict access to #{name}")
        Rollbar.info("Attempting to restrict access to #{name}")
        update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
        update_request.update_revision('patient', name, restricted: true)
        magma_client.update(update_request)
      end
    end

    def unrestrict!(patient)
      name = patient['name']
      logger.warn("Attempting to unrestrict access to #{name}")
      Rollbar.info("Attempting to unrestrict access to #{name}")

      # This code path should be --eventually consistent--  That is to say, we should ensure each operation
      # is idempotent (may need to be repeated), and that the patient is not marked restricted until all other
      # related tasks are complete and the state is consistent.
      mvir1_waiver.release_patient_data(name)

      update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
      update_request.update_revision('patient', name, restricted: false)
      magma_client.update(update_request)
    end

    def mvir1_waiver
      require_relative './metis/mvir1_waiver'
      @mvir1_waiver ||= Mvir1Waiver.new(metis_client: metis_client)
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
