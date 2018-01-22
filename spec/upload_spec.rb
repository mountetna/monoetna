describe "UploadController" do
  include Rack::Test::Methods

  def app
    OUTER_APP
  end

  context "#authorize" do
    it "should authorize an upload" do
      params = {
        file_name: 'wisdom.txt',
        project_name: 'athena',
        user_email: 'metis@ucsf.edu'
      }

      json_post('upload/authorize', params)

      expect(last_response.body).to eq(signature)
    end
  end

  context "#start" do
    it "should start an upload" do
      json_post('upload/start', { })

      expect(last_response.status).to eq(200)
      upload = Upload.first
      resource = Resource.first

      expect(upload).not_to be_nil
      expect(resource).not_to be_nil
    end
  end
end
