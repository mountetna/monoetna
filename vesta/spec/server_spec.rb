describe StatsController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  describe '#global_stats' do
    it 'gets global stats' do
      stats = {
        assay_count: 4,
        byte_count: 30,
        file_count: 3,
        sample_count: 5,
        subject_count: 3,
        user_count: 4
      }
      create(:global_stats, stats)
      get("/stats")

      expect(last_response.status).to eq(200)

      expect(json_body.first).to include(stats)
      expect(json_body.first).to have_key(:recorded_at)
    end
  end

  describe '#project_stats' do
    it 'gets project stats' do
      athena_stats = {
       assay_count: 0,
       byte_count: 10,
       clinical_data_count: 0,
       file_count: 1,
       sample_count: 1,
       subject_count: 1,
       user_count: 1,
      }
      labors_stats = {
       assay_count: 4,
       byte_count: 20,
       clinical_data_count: 0,
       file_count: 2,
       sample_count: 4,
       subject_count: 2,
       user_count: 0
      }
      create(:project_stats, name: "athena", **athena_stats)
      create(:project_stats, name: "labors", **labors_stats)
      get("/stats/projects")

      expect(json_body).to contain_exactly(
       include(
         :assay_count=>0,
         :byte_count=>10,
         :clinical_data_count=>0,
         :file_count=>1,
         :name=>"athena",
         :sample_count=>1,
         :subject_count=>1,
         :user_count=>1),
        include(:assay_count=>4,
         :byte_count=>20,
         :clinical_data_count=>0,
         :file_count=>2,
         :name=>"labors",
         :sample_count=>4,
         :subject_count=>2,
         :user_count=>0)
      )
    end
  end
end

describe ProjectController do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end
  describe '#projects' do
    it 'gets project info' do
      athena_info = {
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
      }
      labors_info = {
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
      }
      create(:project, start_date: Time.now, **athena_info)
      create(:project, start_date: Time.now, **labors_info)
      get("/projects")

      expect(json_body).to contain_exactly(
        include(athena_info),
        include(labors_info)
      )
    end
  end
end
