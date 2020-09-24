require 'date'
require 'logger'
require 'rollbar'
require 'sequel'
require_relative 'helpers'


class Polyphemus
  class Help < Etna::Command
    usage 'List this help'

    def execute
      puts 'Commands:'
      Polyphemus.instance.commands.each do |name, cmd|
        puts cmd.usage
      end
    end
  end

  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'

    def execute(version=nil)
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
    usage 'link_comet_bulk_rna [environment]'

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
        super(**{bucket_name: 'data', model_name: 'rna_seq'}.update(opts))
      end

      def matching_expressions
        [
            [/bulk_RNASeq\/raw\/(?<sample>[^\/]*)\/.*$/, 'raw_fastqs'],
            [/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>_gene_tpms.tsv$/, 'gene_expression'],
            [/bulk_RNASeq\/processed\/(?<sample>[^\/]*)\/\k<sample>_transcript_tpms.tsv$/, 'isoform_expression'],
        ]
      end

      def attribute_options
        { 'raw_fastqs' => { file_collection: true }}
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

  class ApplyMvir1RnaSeqAttributes < Etna::Command
    include WithEtnaClients
    include WithLogger

    attr_reader :environment
    usage 'apply_mvir1rna_seq_attributes [environment]'

    def execute(env = Polyphemus.instance.environment)
      @environment = env
      models = Etna::Clients::Magma::Models.new
      rna_seq = models.build_model('rna_seq')
      attributes = rna_seq.build_template.build_attributes
      add_file_attribute(attributes, 'gene_expression')
      add_file_attribute(attributes, 'isoform_expression')
      add_file_attribute(attributes, 'genome_alignment')
      add_file_attribute(attributes, 'genome_alignment_idx')
      add_file_attribute(attributes, 'transcriptome_alignment')
      add_file_attribute(attributes, 'fusion_gene_junctions')
      add_file_attribute(attributes, 'non_host_reads')
      add_file_attribute(attributes, 'rrna_alignment')
      add_file_attribute(attributes, 'rrna_alignment_idx')
      add_file_attribute(attributes, 'genome_alignment_qc_pdf')
      add_file_attribute(attributes, 'genome_alignment_flagstat')
      add_file_attribute(attributes, 'rrna_alignment_qc_pdf')
      add_file_attribute(attributes, 'rrna_alignment_flagstat')
      add_file_attribute(attributes, 'adapter_trimming_metrics')
      add_file_attribute(attributes, 'rrna_alignment_idx')
      add_file_attribute(attributes, 'rrna_alignment_idx')
      add_file_attribute(attributes, 'rrna_alignment_idx')
      add_string_attribute(attributes, 'raw_fastqs')

      Etna::Clients::Magma::ShallowCopyModelWorkflow.new(
          model_name: 'rna_seq',
          target_project: 'mvir1',
          target_client: magma_client,
          source_models: models
      ).ensure_model_tree('rna_seq')
    end

    def add_file_attribute(attributes, name)
      attributes.build_attribute(name).tap do |attribute|
        attribute.attribute_name = name
        attribute.name = attribute.attribute_name
        attribute.display_name = name.split('_').map(&:capitalize).join(' ')
        attribute.attribute_type = Etna::Clients::Magma::AttributeType::FILE
      end
    end

    def add_string_attribute(attributes, name)
      attributes.build_attribute(name).tap do |attribute|
        attribute.attribute_name = name
        attribute.name = attribute.attribute_name
        attribute.display_name = name.split('_').map(&:capitalize).join(' ')
        attribute.attribute_type = Etna::Clients::Magma::AttributeType::STRING
      end
    end
  end

  module CopyCommand
    def self.included(mod)
      mod.instance_eval do
        include WithEtnaClientsByEnvironment
        include WithLogger

        usage "#{name.snake_case.split(/::/).last} <source_env> <target_env> <source_project> <target_project> <model_name>"

        def execute(source_env, target_env, source_project, target_project, model_name)
          workflow = WORKFLOW.from_api_source(
              model_name: model_name,
              source_project: source_project,
              source_client: environment(source_env).magma_client,
              target_project: target_project,
              target_client: environment(target_env).magma_client,
          )

          workflow.ensure_model_tree(model_name)
        end

        def setup(config)
          super
          Polyphemus.instance.setup_logger
        end
      end
    end
  end

  class ApiCopyModelShallow < Etna::Command
    WORKFLOW = Etna::Clients::Magma::ShallowCopyModelWorkflow
    include CopyCommand
  end

  class ApiCopyModelDeep < Etna::Command
    WORKFLOW = Etna::Clients::Magma::ModelSynchronizationWorkflow
    include CopyCommand
  end

  class ApiAddProject < Etna::Command
    include WithEtnaClientsByEnvironment

    usage 'api_add_project <environment> <project_name>'

    def execute(env, project_name)
      client = environment(env).magma_client
      client.update_model(Etna::Clients::Magma::UpdateModelRequest.new(
          project_name: project_name,
          actions: [Etna::Clients::Magma::AddProjectAction.new]))
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

          mvir1_waiver.restrict_pool_data(pool)

          update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
          update_request.update_revision("#{base_model}_pool", pool, restricted: true)
          magma_client.update(update_request)
        else
          logger.info "#{base_model}_pool #{pool} does not include a restricted patient, relaxing."

          mvir1_waiver.release_pool_data(pool)

          update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
          update_request.update_revision("#{base_model}_pool", pool, restricted: false)
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

      # This code path should be --eventually consistent--  That is to say, we should ensure each operation
      # is idempotent (may need to be repeated), and that the patient is not marked restricted until all other
      # related tasks are complete and the state is consistent.
      mvir1_waiver.release_patient_data(name)

      if patient['restricted']
        Rollbar.info("Attempting to unrestrict access to #{name}")
        update_request = Etna::Clients::Magma::UpdateRequest.new(project_name: project)
        update_request.update_revision('patient', name, restricted: false)
        magma_client.update(update_request)
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

  class FindIpiBulkRnaSeqFolders < Etna::Command
    include WithEtnaClients

    usage "Fetch a list of Metis BulkRNASeq folders for IPI, from integral_data"

    def project
      :mvir1
    end

    def execute
      folders = metis_client.find(
          Etna::Clients::Metis::FindRequest.new(
            project_name: 'ipi',
            bucket_name: 'integral_data',
            params: [Etna::Clients::Metis::FindParam.new(
              attribute: 'name',
              predicate: 'glob',
              value: 'BulkRNASeq/IPI*',
              type: 'folder'
            )])).folders.all
      puts "Found a total of #{folders.length} BulkRNASeq folders in the IPI integral_data bucket"
    end

    def setup(config)
      super
    end
  end

  class IpiAddFlowModel < Etna::Command
    include WithEtnaClients
    include WithLogger

    usage "Add the Flow model to IPI, and migrate data from Sample model to Flow"

    def project_name
      :ipi
    end

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)
    end

    def execute
      require_relative './ipi/migrations/01_flow_model_migration'
      migration = IpiAddFlowModelMigration.new(magma_client: magma_client, magma_crud: magma_crud)
      migration.execute
    end

    def setup(config)
      super
    end
  end

  class LinkIpiFlowWsp < Etna::Command
    include WithEtnaClients
    include WithLogger
    usage 'link_ipi_flow_wsp [environment]'

    attr_reader :environment

    def project_name
      :ipi
    end

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)
    end

    def execute(env = Polyphemus.instance.environment)
      require_relative './ipi/flow_wsp_file_linker'
      @environment = env

      linker = IpiFlowWspLinker.new(magma_crud: magma_crud, metis_client: metis_client, project_name: project_name)
      linker.link_files
    end
  end

  class LinkIpiFlowFcs < Etna::Command
    include WithEtnaClients
    include WithLogger
    usage 'link_ipi_flow_fcs [environment]'

    attr_reader :environment

    def project_name
      :ipi
    end

    def magma_crud
      @magma_crud ||= Etna::Clients::Magma::MagmaCrudWorkflow.new(magma_client: magma_client, project_name: project_name)
    end

    def execute(env = Polyphemus.instance.environment)
      require_relative './ipi/flow_fcs_file_linker'
      @environment = env

      linker = IpiFlowFcsLinker.new(magma_crud: magma_crud, metis_client: metis_client, project_name: project_name)
      linker.link_files
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
