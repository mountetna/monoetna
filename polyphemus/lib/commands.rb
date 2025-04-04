require 'date'
require 'logger'
require 'rollbar'
require 'sequel'
require 'tempfile'
require 'active_support/all'
require_relative 'helpers'
require_relative 'data_processing/xml_dsl'
require_relative 'data_processing/magma_dsl'
require_relative 'data_processing/flow_jo_dsl'
require_relative 'data_eng/argo_workflow_manager'
require_relative 'ipi/process_rna_seq_output'

class Polyphemus
  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'
    string_flags << '--version'

    def execute(version: nil)
      Sequel.extension(:migration)
      db = Polyphemus.instance.db

      if version
        puts "Migrating to version #{version}"
        Sequel::Migrator.run(db, 'db/migrations', target: version.to_i)
      else
        puts 'Migrating to latest'
        Sequel::Migrator.run(db, 'db/migrations')
      end
    end

    def setup(config)
      super
      Polyphemus.instance.setup_db
    end
  end

  class LinkCometBulkRna < Etna::Command
    include WithLogger
    include WithEtnaClients

    attr_reader :environment

    def project_name
      mvir1
    end

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)
    end

    def execute(env = Polyphemus.instance.environment)
      @environment = env

      linker = CometBulkRnaLinker.new(magma_crud: magma_crud, metis_client: metis_client, project_name: project_name)
      linker.link_files
    end

    # Use this once we build out the other attributes / ipi bulk loader.
    # class BulkRnaLinker < Etna::Clients::Magma::FileLinkingWorkflow
    #   def matching_expressions
    #     [
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>\.rsem.genes.results$/, 'gene_expression'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>\.rsem.isoforms.results$/, 'isoform_expression'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>\.trimmed.non_rrna.star.Aligned.sortedByCoout.deduplicated.cram$/, 'genome_alignment'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.non_rrna.star.Aligned.sortedByCoout.deduplicated.cram.crai$/, 'genome_alignment_idx'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.non_rrna.star.Aligned.toTranscriptout.cram$/, 'transcriptome_alignment'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.non_rrna.star.Chimeout.junction$/, 'fusion_gene_junctions'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.non_rrna.star.Unmapout.mate1.fastqoutput\/\k<sample>\/\k<sample>.trimmed.non_rrna.star.Unmapout.mate2.fastq.gz$/, 'non_host_reads'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.rrna.sorted.cram$/, 'rrna_alignment'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.rrna.sorted.cram.crai$/, 'rrna_alignment_idx'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.non_rrna.star.Aligned.sortedByCoout.deduplicated_bamqc.pdf$/, 'genome_alignment_qc_pdf'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.non_rrna.star.Aligned.sortedByCoout.deduplicated.flagstat$/, 'genome_alignment_flagstat'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.rrna.sorted_bamqc.pdf$/, 'rrna_alignment_qc_pdf'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>.trimmed.rrna.sorted.flagstat$/, 'rrna_alignment_flagstat'],
    #         [/\/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/fastp.html$/, 'adapter_trimming_metrics'],
    #     ]
    #   end
    # end

    class CometBulkRnaLinker < Etna::Clients::Magma::FileLinkingWorkflow
      def initialize(**opts)
        super(**{ bucket_name: 'data', model_name: 'rna_seq' }.update(opts))
      end

      def matching_expressions
        [
          [/bulk_RNASeq\/raw\/(?<sample>[^\/]*)\/.*$/, 'raw_fastqs'],
          [/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>_gene_tpms.tsv$/, 'gene_expression'],
          [/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>_transcript_tpms.tsv$/, 'isoform_expression'],
        ]
      end

      def attribute_options
        { 'raw_fastqs' => { file_collection: true } }
      end

      # Subclasses should override this to implement custom logic for how regex matches should match to linking.
      def matches_to_record_identifiers(match_data)
        super(match_data).tap do |identifiers|
          identifiers.update(patient_timepoint_from(match_data['sample']))
          identifiers['rna_seq'] = match_data['sample']
        end
      end

      def revision_for(id, attribute_name, file_path, match_map, record_identifiers)
        super(id, attribute_name, file_path, match_map, record_identifiers).tap do |revision|
          revision['tissue_type'] = 'ETA'
          revision['sequencer'] = 'NovaSeq'
          revision['NovaSeq'] = 'GRCh38_U13369.1_ERCC'
          revision['pipeline_description'] = 'langelier_group'
        end
      end
    end
  end

  class CopyMetisFilesCommand < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger

    def execute(source_env, target_env, source_project, target_project, source_bucket, target_bucket, file_glob_match)
      download_workflow = Etna::Clients::Metis::MetisDownloadWorkflow.new(
        metis_client: environment(source_env).metis_client,
        project_name: source_project,
        bucket_name: source_bucket,
      )

      upload_workflow = Etna::Clients::Metis::MetisUploadWorkflow.new(
        metis_client: environment(target_env).metis_client,
        project_name: target_project,
        bucket_name: target_bucket,
      )

      files = download_workflow.metis_client.find(Etna::Clients::Metis::FindRequest.new(
        project_name: source_project,
        bucket_name: source_bucket,
        params: [
          Etna::Clients::Metis::FindParam.new(
            attribute: 'name',
            predicate: 'glob',
            value: file_glob_match,
          )
        ]
      )).files.all

      logger.info("Found #{files.length} matches for match '#{file_glob_match}'")
      tmpfile = Tempfile.new("download-buffer")
      files.each do |file|
        logger.info("Beginning download of #{file.file_path}")
        download_workflow.do_download(tmpfile.path, file) do |progress|
          # TODO: Maybe add a progress bar?
        end

        logger.info("Beginning upload of #{file.file_path}")
        upload_workflow.do_upload(tmpfile.path, file.file_path) do |progress|
          case progress[0]
          when :error
            logger.warn("Error while uploading: #{progress[1].to_s}")
          else
          end
        end
      end
    end

    def setup(config)
      super
      Polyphemus.instance.setup_logger
    end
  end

  class CopyMetisTreeCommand < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger

    def execute(source_env, target_env, source_project, target_project, source_bucket, target_bucket, folder_name_glob)
      source_metis_client = environment(source_env).metis_client
      dest_metis_client = environment(target_env).metis_client

      folders = source_metis_client.find(Etna::Clients::Metis::FindRequest.new(
        project_name: source_project,
        bucket_name: source_bucket,
        params: [
          Etna::Clients::Metis::FindParam.new(
            type: 'folder',
            attribute: 'name',
            predicate: 'glob',
            value: folder_name_glob,
          )
        ]
      )).folders.all

      logger.info("Found #{folders.length} folders")
      folders.each do |folder|
        dest_metis_client.create_folder(Etna::Clients::Metis::CreateFolderRequest.new(
          project_name: target_project,
          bucket_name: target_bucket,
          folder_path: folder.folder_path
        )
        )
      end
      logger.info("Done")
    end

    def setup(config)
      super
      Polyphemus.instance.setup_logger
    end
  end

  class CopyMagmaRecords < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger

    def execute(source_env, target_env, project_name, *models)
      workflow = Etna::Clients::Magma::RecordSynchronizationWorkflow.new(
        target_client: environment(target_env).magma_client,
        source_client: environment(source_env).magma_client,
        project_name: project_name,
        ignore_update_errors: true,
      )

      if models.empty?
        models = workflow.target_models.model_keys
      end

      logger.info("Copying records from #{source_env} #{project_name} -> #{target_env} #{project_name}")

      models.each do |model|
        workflow.copy_model(model) do |model_name, documents|
          logger.info("Copying #{documents.document_keys.length} #{model_name} records")
        end
      end
    end

    def setup(config)
      super
      Polyphemus.instance.setup_logger
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
      request.attribute_names = ['name', 'consent', 'restricted']
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
          # Magma does not support delete right now, so we restrict the models,
          # but we can unlink and delete the files in metis
          restrict!(patient, delete_metis_files: true)
        elsif should_be_restricted
          restrict!(patient)
        elsif !should_be_restricted
          unrestrict!(patient)
        end
      end

      # Check all pools and verify status
      cascade_to_cytof_pools
      cascade_to_sc_rna_pools
    end

    def cascade_to_cytof_pools
      cascade_to_pool('cytof')
    end

    def cascade_to_sc_rna_pools
      cascade_to_pool('sc_rna_seq')
    end

    def cascade_to_pool(base_model)
      request = Etna::Clients::Magma::QueryRequest.new(project_name: project)
      request.query = [base_model,
        ['timepoint', 'patient', 'restricted', '::true'],
        '::all', "#{base_model}_pool", '::identifier']
      all_restricted_pools = magma_client.query(request).answer.map { |r| r[1] }.sort.uniq

      request.query = ["#{base_model}_pool", '::all', '::identifier']
      all_pools = magma_client.query(request).answer.map { |r| r[1] }.sort.uniq
      all_pools.each do |pool|
        if all_restricted_pools.include? pool
          logger.info "#{base_model}_pool #{pool} includes a restricted patient, restricting."

          begin
            mvir1_waiver.restrict_pool_data(pool)
          rescue Etna::Error => e
            logger.log_error(e)
          end

          update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
          update_request.update_revision("#{base_model}_pool", pool, restricted: true)
          magma_client.update_json(update_request)
        else
          logger.info "#{base_model}_pool #{pool} does not include a restricted patient, relaxing."

          begin
            mvir1_waiver.release_pool_data(pool)
          rescue Etna::Error => e
            logger.log_error(e)
          end

          update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
          update_request.update_revision("#{base_model}_pool", pool, restricted: false)
          magma_client.update_json(update_request)
        end
      end
    end

    def restrict!(patient, delete_metis_files: false)
      name = patient['name']

      begin
        mvir1_waiver.restrict_patient_data(name, delete_metis_files)
      rescue Etna::Error => e
        logger.log_error(e)
      end

      unless patient['restricted']
        logger.warn("Attempting to restrict access to #{name}")
        Rollbar.info("Attempting to restrict access to #{name}")
        update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
        update_request.update_revision('patient', name, restricted: true)
        magma_client.update_json(update_request)
      end
    end

    def unrestrict!(patient)
      name = patient['name']
      logger.warn("Attempting to unrestrict access to #{name}")

      # This code path should be --eventually consistent--  That is to say, we should ensure each operation
      # is idempotent (may need to be repeated), and that the patient is not marked restricted until all other
      # related tasks are complete and the state is consistent.
      begin
        mvir1_waiver.release_patient_data(name)
      rescue Etna::Error => e
        logger.log_error(e)
      end

      if patient['restricted']
        Rollbar.info("Attempting to unrestrict access to #{name}")
        update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
        update_request.update_revision('patient', name, restricted: false)
        magma_client.update_json(update_request)
      end
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

  # this is still useful, maybe convert???
  class LinkIpiFlowWsp < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'link_ipi_flow_wsp <environment>'

    def project_name
      :ipi
    end

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: @environ.magma_client, project_name: project_name)
    end

    def execute(env)
      require_relative './ipi/flow_wsp_file_linker'
      @environ = environment(env)

      linker = IpiFlowWspLinker.new(
        magma_crud: magma_crud,
        metis_client: @environ.metis_client,
        project_name: project_name)
      linker.link_files
    end
  end

  class RelinkFileCollectionAttribute < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
        magma_client: @environ.magma_client, project_name: project_name)
    end

    string_flags << "--env"
    string_flags << "--project_name"
    string_flags << "--model_name"
    string_flags << "--attribute_name"

    def execute(
      env:,
      project_name:,
      model_name:,
      attribute_name:
    )
      magma_client = environment(env).magma_client
      models = magma_models(magma_client, project_name)

      if (attribute = models.model(model_name).template.attributes.attribute(attribute_name))
        unless attribute.attribute_type == Etna::Clients::Magma::AttributeType::FILE_COLLECTION
          unless (migrated_attribute = models.model(model_name).template.attributes.attribute("#{attribute_name}_migrated"))
            magma_client.update_model(
              Etna::Clients::Magma::UpdateModelRequest.new(
                project_name: project_name,
                actions: [
                  Etna::Clients::Magma::AddAttributeAction.new(
                    model_name: model_name,
                    attribute_name: "#{attribute_name}_migrated",
                    type: Etna::Clients::Magma::AttributeType::FILE_COLLECTION,
                    description: attribute.description,
                    display_name: attribute.display_name,
                  )
                ]
              )
            )

            # models = magma_models(magma_client, project_name)
            # migrated_attribute = models.model(model_name).template.attributes.attribute("#{attribute_name}_migrated")
          end
        end

        puts "hiding #{attribute_name}"
        magma_client.update_model(
          Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project_name,
            actions: [Etna::Clients::Magma::UpdateAttributeAction.new(
              model_name: model_name,
              attribute_name: attribute_name,
              hidden: true,
            )]
          )
        )

        puts "renaming #{attribute_name} -> #{attribute_name}_old"
        magma_client.update_model(
          Etna::Clients::Magma::UpdateModelRequest.new(
            project_name: project_name,
            actions: [Etna::Clients::Magma::RenameAttributeAction.new(
              model_name: model_name,
              attribute_name: attribute_name,
              new_attribute_name: "#{attribute_name}_old"
            )],
          )
        )
      end

      puts "renaming #{attribute_name}_migrate -> #{attribute_name}s"
      magma_client.update_model(
        Etna::Clients::Magma::UpdateModelRequest.new(
          project_name: project_name,
          actions: [Etna::Clients::Magma::RenameAttributeAction.new(
            model_name: model_name,
            attribute_name: "#{attribute_name}_migrated",
            new_attribute_name: attribute_name
          )]
        )
      )
    end

    def magma_models(magma_client, project_name)
      magma_client.retrieve(Etna::Clients::Magma::RetrievalRequest.new(
        project_name: project_name,
        model_name: "all",
        attribute_names: "all",
        record_names: [],
      )).models
    end
  end

  # Keep until ipi is 'done' ???
  class LinkIpiFlowFcs < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'link_ipi_flow_fcs <environment>'

    def project_name
      :ipi
    end

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
        magma_client: @environ.magma_client, project_name: project_name)
    end

    def execute(env)
      require_relative './ipi/flow_fcs_file_linker'
      @environ = environment(env)

      linker = IpiFlowFcsLinker.new(
        magma_crud: magma_crud,
        metis_client: @environ.metis_client,
        project_name: project_name)
      linker.link_files
    end
  end

  # First step of workflow
  class IpiCopyFlowToIntegralDataset < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'ipi_copy_flow_to_integral_dataset <environment> <source_bucket> <source_folder>'

    def execute(env, source_bucket, source_folder)
      require_relative './ipi/flow_populate_integral_data'
      metis_client = environment(env).metis_client

      integral_flow = IpiFlowPopulateIntegralData.new(
        metis_client: metis_client,
        source_bucket_name: source_bucket,
        source_folder_name: source_folder
      )
      integral_flow.copy_files
    end
  end

  class SetFileAttributesToBlank < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'set_file_attributes_to_blank <environment> <project_name> <model_names>'

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
        magma_client: @environ.magma_client,
        project_name: @project_name)
    end

    def execute(env, project_name, *model_names)
      @environ = environment(env)
      @project_name = project_name

      blanker = Etna::Clients::Magma::FileAttributesBlankWorkflow.new(
        magma_crud: magma_crud,
        model_names: model_names,
        project_name: project_name)
      blanker.set_file_attrs_blank
    end
  end

  class UpdateAttributesFromCsv < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'update_attributes_from_csv <environment> <project_name> <model_name> <filepath>'

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(
        magma_client: @environ.magma_client,
        project_name: @project_name)
    end

    def execute(env, project_name, model_name, filepath)
      @environ = environment(env)
      @project_name = project_name

      update_attributes_workflow = Etna::Clients::Magma::UpdateAttributesFromCsvWorkflowSingleModel.new(
        magma_crud: magma_crud,
        project_name: project_name,
        model_name: model_name,
        filepath: filepath)
      update_attributes_workflow.update_attributes
    end
  end

  class LoadRnaSeqProcessedData < Etna::Command
    include WithEtnaClientsByEnvironment
    include WithLogger
    usage 'load_rna_seq_processed_data <env> <project_name> <file_path> [--execute]'
    boolean_flags << '--execute'

    def execute(env, project_name, file_path, execute: false)
      @environ = environment(env)
      @project_name = project_name

      process_rna_seq_output = ProcessRnaSeqOutput.new(
        magma_client: @environ.magma_client,
        project_name: project_name,
        file_path: file_path,
        execute: execute)

      process_rna_seq_output.process_rna_seq
    end
  end

  class CopyGNEPoolData < Etna::Command
    include WithEtnaClients
    include WithLogger

    def execute
      workflow = Etna::Clients::Metis::SyncMetisDataWorkflow.new(
        metis_client: metis_client, logger: logger,
        project_name: 'mvir1', bucket_name: 'GNE_redacted_data',
        filesystem: filesystem)

      workflow.copy_directory("data", "/Upload/processed/GNE_redacted_data", "/Upload/processed", nil)
      logger.info("Done")
    end

    def filesystem
      aspera_comet = Polyphemus.instance.config(:aspera_comet)
      @filesystem ||= Etna::Filesystem::GeneAsperaCliFilesystem.new(**aspera_comet)
    end

    def setup(config)
      super
      Polyphemus.instance.setup_logger
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
      Polyphemus.instance.setup_db
      Polyphemus.instance.setup_sequel
      Polyphemus.instance.setup_ssh
    end
  end

  class RunJob < Etna::Command
    include WithEtnaClients

    def setup(config)
      super
      Polyphemus.instance.setup_logger
      Polyphemus.instance.setup_db
      Polyphemus.instance.setup_sequel
    end

    def execute(job_name, config_id, version_number)
      # Retrieve the config from polyphemus
      # We need to use the db here to get the decrypted secrets
      config = Polyphemus::Config.current.where(
        config_id: config_id,
        version_number: version_number
      ).first

      runtime_config = Polyphemus::RuntimeConfig.where(
        config_id: config_id,
      ).first

      @token = janus_client.generate_token(
          'task',
          signed_nonce: nil,
          project_name: config.project_name
      )

      # reset clients so they use task token
      reset_clients!

      # Instantiate the job and run it
      # Dynamically load the job class from the job_name
      job = Kernel.const_get("#{job_name}_job".camelize.to_sym)

      job.new(
        @token,
        config.with_secrets,
        runtime_config.as_json
      ).execute
    end

  end

  class GetRuntimeMetadata < Etna::Command
    include WithEtnaClients

    def setup(config)
      super
      Polyphemus.instance.setup_logger
      Polyphemus.instance.setup_db
      Polyphemus.instance.setup_sequel
    end

    def execute(run_id, workflow_name, workflow_namespace)

      workflow_json = %x{ argo get --output=json #{workflow_name} -n #{workflow_namespace} }
      raw_output = %x{ argo logs #{workflow_name} }

      # We need to fetch the project name
      run = Polyphemus::Run.where(
        run_id: run_id
      ).first

      config = Polyphemus::Config.current.where(
        config_id: run.config_id,
        version_number: run.version_number
      ).first
      
      # Parse the workflow_json and just extract status
      workflow_data = JSON.parse(workflow_json)
      status = workflow_data["status"]
      name = workflow_data.dig("metadata","name")

      # remove control characters and workflow name from output logs
      output = raw_output.gsub(/\e\[\d+m/,'').gsub(/^#{name}[\-\w]*: /,'')

      updates = {
        orchestrator_metadata: status,
        output: output
      }
      polyphemus_client.update_run(config.project_name, run_id, updates)
    end

  end

  class IntervalScheduler < Etna::Command
    include WithLogger
  
    usage 'Continuously polls eligible runtime configs and submits Argo workflows' 

    SLEEP_INTERVAL = 60 * 5 # 5 Minutes

    def setup(config)
      super
      Polyphemus.instance.setup_logger
      Polyphemus.instance.setup_db
      Polyphemus.instance.setup_sequel
    end
  
    def execute
      loop do
        begin
          eligible_runtime_configs = Polyphemus::RuntimeConfig.eligible_runtime_configs
          logger.info("Found #{eligible_runtime_configs.count} eligible runtime configs for scheduling...")
          eligible_runtime_configs.each do |runtime_config|
            config = Polyphemus::Config.current.where(
              config_id: runtime_config.config_id,
            ).first
            Polyphemus::ArgoWorkflowManager.submit_workflow(config)
            sleep SLEEP_INTERVAL
          end
        rescue StandardError => e
          logger.error("Error in scheduler loop: #{e}")
          sleep SLEEP_INTERVAL
        end
      end
    end
  end
end
