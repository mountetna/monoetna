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

      header(
        *Etna::TestAuth.token_header(
          email: 'metis@ucsf.edu', perm: 'e:athena'
        )
      )

      json_post('authorize/upload', params)
      
      url = last_response.body
      uri = URI(url)

      hmac_params = Rack::Utils.parse_nested_query(uri.query)

      expect(last_response.status).to eq(200)
      expect(uri.path).to eq('/upload')
      expect(hmac_params['X-Etna-Id']).to eq('metis')
      expect(hmac_params['X-Etna-Project-Name']).to eq(params[:project_name])
      expect(hmac_params['X-Etna-File-Name']).to eq(params[:file_name])

      # now post it
      url = url.sub(%r!https://example.org/!,'')
      header('Authorization', nil)
      
      json_post(url, action: 'start')

      #expect(last_response.status).to eq(200)
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
