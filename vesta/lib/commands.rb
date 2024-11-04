require_relative 'helpers'

class Vesta
  class Migrate < Etna::Command
    include WithLogger

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
    end
  end

  class CollectDLStats < Etna::Command
    include WithEtnaClients
    include WithLogger

    usage 'Collect Data Library-wide usage stats'

    def execute
      # get all projects
      janus_stats = janus_client.get_project_stats

      # byte and file counts
      file_count_by_project = metis_client.get_file_count_by_project
      byte_count_by_project = metis_client.get_byte_count_by_project

      project_stats = janus_stats[:projects].map do |proj|
        # skip resource projects
        # ToDo: Make it not skip community projects
        next if proj[:resource] == true and proj[:requires_agreement] == false

        proj_name = proj[:project_name]

        # get project template
        begin
          puts "Collecting Magma counts for #{proj_name}"
          models = magma_client
            .retrieve(project_name: proj_name, model_name: 'all', attribute_names: [], record_names: [])
            .models
            .model_keys
            .to_set

          subject_count = query_model_counts(proj_name, @data_types[:subject_models].to_set, models)
          sample_count = query_model_counts(proj_name, @data_types[:sample_models].to_set, models)
          assay_count = query_model_counts(proj_name, @data_types[:assay_models].to_set, models)
          clinical_data_count = query_model_counts(proj_name, @data_types[:clinical_models].to_set, models)
        rescue => e
          puts "Error collecting Magma counts for #{proj_name}: #{e}"
          next
        end

        {
          name: proj_name,
          user_count: proj[:user_count],
          file_count: file_count_by_project[proj_name.to_sym] || 0,
          byte_count: byte_count_by_project[proj_name.to_sym] || 0,
          subject_count: subject_count,
          sample_count: sample_count,
          assay_count: assay_count,
          clinical_data_count: clinical_data_count,
        }
      end.compact

      puts "DONE collecting project counts"

      global_stats = Hash.new(0)
      project_stats.each do |item|
        item.each do |k, v|
          next if [:name, :user_count, :clinical_data_count].include?(k)
          global_stats[k] += v
        end
      end

      # top-level user count dedupes users belonging to multiple projects
      global_stats[:user_count] = janus_stats[:user_count]
      puts "DONE generating global stats"

      Vesta.instance.db.transaction do
        # save projects stats
        Vesta::ProjectStats.multi_insert(project_stats)
        # save global stats
        Vesta::GlobalStats.create(global_stats)
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

      count || 0
    end

    def setup(config)
      super
      @data_types = Vesta.instance.config(:data_types)
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
        proj_name = proj[:project_name]

        # skip resource projects, ToDo: don't skip community projects
        if proj[:resource] == true and proj[:requires_agreement] == false
          puts "Skipping collecting project info for #{proj_name}"
          next
        end

        puts "Start collecting project info for #{proj_name}"

        begin
          project_info = retrieve_project_info(proj_name)
          project_model_names = retrieve_model_names(proj_name)

          pi_profiles = proj[:principal_investigators].map do |janus_pi|
            begin
              retrieve_ucsf_profile(janus_pi)
            rescue => e
              puts "Error retreiving UCSF profile for #{janus_pi[:name]}: #{e}"
              {
                name: janus_pi[:name],
                email: janus_pi[:email],
                profile_url: nil,
                title: nil,
                photo_url: nil,
              }
            end
          end

          # insert project data one-by-one (vs multi_insert)
          # in case any one project has null vals for non-nullable attrs
          Vesta::Project.update_or_create(
            {
              name: proj_name,
            },
            full_name: proj[:project_name_full],
            description: project_info[:description],
            funding_source: project_info[:funding_source],
            principal_investigators: pi_profiles,
            status: project_info[:project_status],
            type: project_info[:project_type],
            species: project_info[:species],
            start_date: parse_start_date(project_info[:start_date]),
            theme: project_info[:theme],
            data_collection_complete: project_info[:completed] != nil ? project_info[:completed] : false,
            data_types: project_model_names,
          )
        rescue => e
          puts "Error collecting project info for #{proj_name}: #{e}"
          next
        end

        puts "Successfully collected project info for #{proj_name}"
      end
    end

    PROJECT_ATTRIBUTES=[
        'completed',
        'description',
        'funding_source',
        'project_status',
        'project_type',
        'species',
        'start_date',
        'theme'
    ]

    def retrieve_project_info(project_name)
      response = magma_client.retrieve(
        project_name: project_name,
        model_name: 'project',
        attribute_names: PROJECT_ATTRIBUTES,
        record_names: 'all',
      )
      # There should only ever be 1 project record, but its id is not enforced
      docs = response.models.raw["project"]["documents"]
      data = docs[docs.keys().first]
      data.transform_keys(&:to_sym)
    end

    def retrieve_model_names(project_name)
      response = magma_client.retrieve(
        project_name: project_name,
        model_name: 'all',
        attribute_names: [],
        record_names: [],
      )

      response.models.raw.keys.map(&:to_sym)
    end

    def retrieve_ucsf_profile(janus_pi)
      uri = URI('https://api.profiles.ucsf.edu/json/v2/')
      params = {
        ProfilesURLName: janus_pi[:email].split('@')[0],
        source: 'datalibrary.ucsf.edu',
      }
      uri.query = URI.encode_www_form(params)

      res = Net::HTTP.get_response(uri)

      return janus_pi.merge({
        profile_url: nil,
        title: nil,
        photo_url: nil,
      }) unless res.is_a?(Net::HTTPSuccess)

      json = JSON.parse(res.body, symbolize_names: true)[:Profiles][0]

      return janus_pi.merge({
        profile_url: json[:ProfilesURL],
        title: json[:Title],
        photo_url: json[:PhotoURL],
      })
    end

    def parse_start_date(start_date)
      start_date ? Date.strptime(start_date, '%Y') : nil
    end

    def setup(config)
      super
      @data_types = Vesta.instance.config(:data_types)
      Vesta.instance.setup_db
      Vesta.instance.setup_sequel
      Vesta.instance.setup_logger
    end
  end
end
