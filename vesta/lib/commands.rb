require 'date'
require 'logger'
require 'rollbar'
require 'sequel'
require 'tempfile'
require 'active_support/all'
require 'etna/command'
require_relative 'helpers'


class Vesta
  class Migrate < Etna::Command
    usage 'Run migrations for the current environment.'
    string_flags << '--version'

    def execute(version: nil)
      Sequel.extension(:migration)
      db = Vesta.instance.db

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
      Vesta.instance.setup_db
    end
  end

  class Console < Etna::Command
    usage 'Open a console with a connected Vesta instance.'

    def execute
      require 'irb'
      ARGV.clear
      IRB.start
    end

    def setup(config)
      super
      Vesta.instance.setup_db
      Vesta.instance.setup_sequel
      Vesta.instance.setup_ssh
    end
  end

  class CollectDLStats < Etna::Command
    include WithEtnaClients
    include WithLogger

    usage 'Collect Data Library-wide usage stats'

    def execute
      # get all projects
      janus_stats_by_project = janus_client.get_project_stats

      # byte and file counts
      file_count_by_project = metis_client.get_file_count_by_project
      byte_count_by_project = metis_client.get_byte_count_by_project

      project_stats = janus_stats_by_project[:projects].map do |proj|
        # skip resource projects
        next if proj[:resource] == true

        proj_name = proj[:project_name]

        # get project template
        begin
          models = magma_client
            .retrieve(project_name: proj_name, model_name: 'all', attribute_names: [], record_names: [])
            .models
            .model_keys
            .to_set

          subject_count = query_model_counts(proj_name, @magma_config[:subject_models].to_set, models)
          sample_count = query_model_counts(proj_name, @magma_config[:sample_models].to_set, models)
          assay_count = query_model_counts(proj_name, @magma_config[:assay_models].to_set, models)
        rescue => e
          # log error
          next
        end

        {
          name: proj_name,
          user_count: proj[:user_count],
          file_count: file_count_by_project[proj_name] || 0,
          byte_count: byte_count_by_project[proj_name] || 0,
          subject_count: subject_count,
          sample_count: sample_count,
          assay_count: assay_count,
        }
      end

      Vesta.instance.db.transaction do
        # save projects stats
        Vesta::ProjectStats.multi_insert(project_stats.compact)
        # save global stats
        Vesta::GlobalStats.create(user_count: janus_stats_by_project[:user_count])
      end
    end

    def query_model_counts(project_name, allowed_models, existing_models)
      count = nil

      allowed_models.intersection(existing_models).each do |model|
        query_request = Etna::Clients::Magma::QueryRequest.new(
          project_name: project_name,
          query: [model, "::count"],
        )

        response_count = magma_client.query(query_request).answer

        count = 0 if count.nil?
        count += response_count
      end

      count
    end

    def setup(config)
      super
      @magma_config = Vesta.instance.config(:magma)
      Vesta.instance.setup_db
      Vesta.instance.setup_sequel
      Vesta.instance.setup_logger
    end
  end

  class CollectProjectInfo < Etna::Command
    include WithEtnaClients
    include WithLogger

    usage 'Collect standard project info for all projects'

    def execute
      # get all projects
      janus_stats_by_project = janus_client.get_project_stats

      janus_stats_by_project[:projects].each do |proj|
        # skip resource projects
        next if proj[:resource] == true

        proj_name = proj[:project_name]

        begin
          project_info = retrieve_project_info(proj_name)

          # insert project data one-by-one (vs multi_insert)
          # in case any one project has null vals for non-nullable attrs
          Vesta::Project.update_or_create(
            {
              name: project_info[:name],
            },
            full_name: proj[:project_name_full],
            description: project_info[:description],
            funding_source: project_info[:funding_source],
            principal_investigator: project_info[:principal_investigator],
            status: project_info[:project_status],
            type: project_info[:project_type],
            species: project_info[:species],
            start_date: project_info[:start_date],
            theme: project_info[:theme],
            data_collection_complete: project_info[:completed],
          )
        rescue => e
          # log error
          next
        end
      end
    end

    def retrieve_project_info(project_name)
      response = magma_client.retrieve(
        project_name: project_name,
        model_name: 'project',
        attribute_names: [
          'completed', 'description', 'funding_source',
          'principal_investigator', 'project_status', 'project_type',
          'species', 'start_date', 'theme',
        ],
        record_names: 'all',
      )

      data = response.models.raw["project"]["documents"][project_name]
      data.transform_keys(&:to_sym)
    end

    def setup(config)
      super
      @magma_config = Vesta.instance.config(:magma)
      Vesta.instance.setup_db
      Vesta.instance.setup_sequel
      Vesta.instance.setup_logger
    end
  end
end