require_relative '../lib/commands'

describe 'Vesta Commands' do
  projects = {
    user_count: 4,
    projects: [
      {
        project_name: "athena",
        project_name_full: "Athena",
        resource: false,
        user_count: 1,
        principal_investigators: [ {
          name: "Metis",
          email: "metis@headofzeus.org"
        } ]
      },
      {
        project_name: "labors",
        project_name_full: "Labors",
        resource: false,
        user_count: 0,
        principal_investigators: [
          {
            name: "Eurystheus",
            email: "eurystheus@twelve-labors.org"
          },
          {
            name: "Hera",
            email: "hera@olympus.org"
          }
        ]
      }
    ]
  }

  describe Vesta::CollectDLStats do
    subject(:collect_global_stats) {
      described_class.new
    }

    it "collects data from other apps" do
      file_counts = {
        athena: 1,
        labors: 2,
      }

      bytes = {
        athena: 10,
        labors: 20
      }

      stub_options("https://metis.test/")

      stub_json("https://janus.test/api/stats/projects", projects)
      stub_json("https://metis.test/api/stats/files", file_counts)
      stub_json("https://metis.test/api/stats/bytes", bytes)

      stub_models("athena", [ "olympian", "blood" ])
      stub_query_count("athena", "olympian", 1)
      stub_query_count("athena", "blood", 1)

      stub_models("labors", [ "victim", "form", "census" ])
      stub_query_count("labors", "victim", 2)
      stub_query_count("labors", "form", 4)
      stub_query_count("labors", "census", 4)

      command = Vesta::CollectDLStats.new

      command.setup(Vesta.instance.instance_variable_get("@config"))

      command.execute

      expect(Vesta::GlobalStats.count).to eq(1)
      expect(Vesta::ProjectStats.count).to eq(2)
      expect(Vesta::GlobalStats.first.to_hash).to include(
        assay_count: 4,
        byte_count: 30,
        file_count: 3,
        sample_count: 5,
        subject_count: 3,
        user_count: 4,
      )
      expect(Vesta::ProjectStats[name: "athena"].to_hash).to include(
       assay_count: 0,
       byte_count: 10,
       clinical_data_count: 0,
       file_count: 1,
       sample_count: 1,
       subject_count: 1,
       user_count: 1,
      )
      expect(Vesta::ProjectStats[name: "labors"].to_hash).to include(
       assay_count: 4,
       byte_count: 20,
       clinical_data_count: 0,
       file_count: 2,
       sample_count: 4,
       subject_count: 2,
       user_count: 0,
      )
    end
  end

  describe Vesta::CollectProjectInfo do
    subject(:collect_global_stats) {
      described_class.new
    }

    def stub_retrieve_project_info(project, response)
      stub_retrieve(
        {
          "project_name":project,
          "model_name":"project",
          "attribute_names": Vesta::CollectProjectInfo::PROJECT_ATTRIBUTES,
          "record_names":"all"
        },
        {
          models: {
            project: {
              documents: { project => response }
            }
          }
        }
      )
    end

    def stub_profile(name, title)
      stub_request(
        :get,
        "https://api.profiles.ucsf.edu/json/v2/?ProfilesURLName=#{name}&source=datalibrary.ucsf.edu"
      ).to_return(
        status: 200,
        body: {
          Profiles: [
            {
              ProfilesURL: "https://profiles.something/#{name}",
              PhotoURL: "https://photos.something/#{name}",
              Title: title
            }
          ]
        }.to_json,
        headers: {}
      )
    end


    it "collects project info from janus and magma" do
      stub_json("https://janus.test/api/stats/projects", projects)

      now = Time.now.iso8601

      stub_retrieve_project_info("athena", {
        name: "athena",
        completed: false,
        description: "Birth Athena",
        funding_source: "Zeus",
        project_status: "in capite",
        project_type: "community",
        species: "olympian",
        start_date: now,
        theme: "fetal"
      })
      stub_retrieve_project_info("labors", {
        name: "labors",
        completed: false,
        description: "Ten labors plus two assigned to Hercules as penance.",
        funding_source: "Hera",
        project_status: "in media res",
        project_type: "community",
        species: "lion,hydra",
        start_date: now,
        theme: "taxonomy"
      })
      stub_models("athena", [ "olympian", "blood" ])
      stub_models("labors", [ "victim", "form", "census" ])
      stub_profile("metis", "Titaness")
      stub_profile("eurystheus", "King")
      stub_profile("hera", "Queen of Heaven")

      command = Vesta::CollectProjectInfo.new

      command.setup(Vesta.instance.instance_variable_get("@config"))

      command.execute

      expect(Vesta::Project.count).to eq(2)
      expect(Vesta::Project[name: "athena"].as_json).to include(
       name: "athena",
       data_collection_complete: false,
       data_types: [ "olympian", "blood" ],
       description: "Birth Athena",
       full_name: "Athena",
       funding_source: "Zeus",
       principal_investigators: [{:email=>"metis@headofzeus.org", :name=>"Metis", :photo_url=>"https://photos.something/metis", :profile_url=>"https://profiles.something/metis", :title=>"Titaness"}],
       species: "olympian",
       status: "in capite",
       theme: "fetal",
       type: "community"
      )
      expect(Vesta::Project[name: "labors"].as_json).to include(
       name: "labors",
       data_collection_complete: false,
       data_types: ["victim", "form", "census"],
       description: "Ten labors plus two assigned to Hercules as penance.",
       full_name: "Labors",
       funding_source: "Hera",
       principal_investigators: [
         {email: "eurystheus@twelve-labors.org", name: "Eurystheus", photo_url: "https://photos.something/eurystheus", profile_url: "https://profiles.something/eurystheus", title: "King"},
         {email: "hera@olympus.org", name: "Hera", photo_url: "https://photos.something/hera", profile_url: "https://profiles.something/hera", title: "Queen of Heaven"}
       ],
       species: "lion,hydra",
       status: "in media res",
       theme: "taxonomy",
       type: "community"
      )
    end
  end

  describe Vesta::Console do
    subject(:console) { described_class.new.execute }
    before do
      require 'irb'
      allow(ARGV).to receive(:clear)
      allow(IRB).to receive(:start)
    end

    it 'calls ARGV and IRB' do
      console

      expect(ARGV).to have_received(:clear).once
      expect(IRB).to have_received(:start).once
    end
  end
end
