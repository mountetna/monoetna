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
            Polyphemus::ArgoWorkflowManager.submit_workflow(runtime_config.config)
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
