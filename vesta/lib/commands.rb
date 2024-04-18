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

  class CollectGlobalStats < Etna::Command
    include WithEtnaClients
    include WithLogger

    usage 'Collect Data Library-wide usage stats'

    def execute(magma_subject_models, magma_sample_models, magma_assay_models)
      # TODO: do all business logic here and separate each stats endpoint into separate stats

      # get all projects
      janus_stats_per_project = janus_client.get_stats

      # byte and file counts
      metis_stats_per_project = metis_client.get_stats

      janus_stats_per_project[:projects].each do |proj|
        # skip resource projects
        next if proj[:resource] == true

        proj_name = proj[:project_name]

        subject_count = query_magma_count_for_possible_models(proj_name, magma_subject_models)
        sample_count = query_magma_count_for_possible_models(proj_name, magma_sample_models)
        assay_count = query_magma_count_for_possible_models(proj_name, magma_assay_models)
      end

      # Vesta.instance.
      # save to DB
    end

    def query_magma_count_for_possible_models(project_name, models)
      count = nil

      models.each do |model|
        query_request = Etna::Clients::Magma::QueryRequest.new(
          [model, "::count"],
          project_name,
          "asc",
          1,
          100
        )

        response = magma_client.query(query_request)

        # if response.type == "Array"
        #   count = response.something

        #   break
        # end
      end

      count
    end

    def setup(config)
      super
      Vesta.instance.setup_db
      Vesta.instance.setup_sequel
      Vesta.instance.setup_logger
    end
  end
end
