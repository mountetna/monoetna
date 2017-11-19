describe "UploadController" do
  include Rack::Test::Methods

  def app
    OUTER_APP
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
